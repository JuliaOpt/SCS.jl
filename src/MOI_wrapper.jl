using MathOptInterface
const MOI = MathOptInterface
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const MOIU = MOI.Utilities

include("matrix.jl")
include("geometric.jl")

"""
    ADMMIterations()

The number of ADMM iterations completed during the solve.
"""
struct ADMMIterations <: MOI.AbstractModelAttribute
end

MOI.is_set_by_optimize(::ADMMIterations) = true

mutable struct MOISolution
    ret_val::Int
    raw_status::String
    objective_value::Float64
    dual_objective_value::Float64
    objective_constant::Float64
    solve_time_sec::Float64
    iterations::Int
end
MOISolution() = MOISolution(0, # SCS_UNFINISHED
                            "", NaN, NaN, NaN,
                            0.0, 0)

function _managed_matrix(A::SparseMatrixCSRtoCSC{T}) where T
    final_touch(A)
    return ManagedSCSMatrix{T}(A.m, A.n, A.nzval, A.rowval, A.colptr)
end

# This is tied to SCS's internal representation
struct ConeData
    qa::Vector{Int} # array of second-order cone constraints
    sa::Vector{Int} # array of semi-definite constraints
    p::Vector{Float64} # array of power cone params
    ConeData() = new(Int[], Int[], Float64[])
end


const CONE_TYPES = (MOI.Zeros, MOI.Nonnegatives, MOI.SecondOrderCone, MOI.PositiveSemidefiniteConeTriangle, MOI.ExponentialCone, MOI.DualExponentialCone, MOI.PowerCone{Float64}, MOI.DualPowerCone{Float64})
const Form{T} = GeometricConicForm{
    Float64,
    SparseMatrixCSRtoCSC{T},
    typeof(CONE_TYPES)
}

mutable struct Optimizer <: MOI.AbstractOptimizer
    cone::ConeData
    data::Union{Form{Int}, Form{Int32}}
    sol::MOISolution
    silent::Bool
    options::Dict{Symbol, Any}
    function Optimizer(; kwargs...)
        optimizer = new(
            ConeData(),
            GeometricConicForm{Float64, SparseMatrixCSRtoCSC{Int}}(CONE_TYPES),
            MOISolution(), false, Dict{Symbol, Any}()
        )
        for (key, value) in kwargs
            MOI.set(optimizer, MOI.RawParameter(String(key)), value)
        end
        return optimizer
    end

end

MOI.get(::Optimizer, ::MOI.SolverName) = "SCS"

function MOI.set(optimizer::Optimizer, param::MOI.RawParameter, value)
    # TODO(odow): remove warning in future version.
    if !(param.name isa String)
        Base.depwarn(
            "passing `$(param.name)` to `MOI.RawParameter` as type " *
            "`$(typeof(param.name))` is deprecated. Use a string instead.",
            Symbol("MOI.set")
        )
    end
    optimizer.options[Symbol(param.name)] = value
end
function MOI.get(optimizer::Optimizer, param::MOI.RawParameter)
    # TODO(odow): remove warning in future version.
    if !(param.name isa String)
        Base.depwarn(
            "passing $(param.name) to `MOI.RawParameter` as type " *
            "$(typeof(param.name)) is deprecated. Use a string instead.",
            Symbol("MOI.get")
        )
    end
    return optimizer.options[Symbol(param.name)]
end

MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool)
    optimizer.silent = value
end
MOI.get(optimizer::Optimizer, ::MOI.Silent) = optimizer.silent

MOI.is_empty(optimizer::Optimizer) = MOI.is_empty(optimizer.data)
function MOI.empty!(optimizer::Optimizer)
    empty!(optimizer.cone.qa)
    empty!(optimizer.cone.sa)
    empty!(optimizer.cone.p)
    MOI.empty!(optimizer.data)
    optimizer.sol.ret_val = 0
end

function MOI.supports(::Optimizer,
                      ::Union{MOI.ObjectiveSense,
                              MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}})
    return true
end

function MOI.supports(::Optimizer, ::MOI.VariablePrimalStart,
                      ::Type{MOI.VariableIndex})
    return true
end

function MOI.supports(::Optimizer,
                      ::Union{MOI.ConstraintPrimalStart,
                              MOI.ConstraintDualStart},
                      ::Type{<:MOI.ConstraintIndex})
    return true
end

function MOI.supports_constraint(optimizer::Optimizer, F::Type{MOI.VectorAffineFunction{Float64}}, S::Type{<:MOI.AbstractVectorSet})
    return MOI.supports_constraint(optimizer.data, F, S)
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kws...)
    linear_solver, _ = sanitize_SCS_options(dest.options)
    T = scsint_t(linear_solver)
    if !(dest.data isa Form{T})
        dest.data = GeometricConicForm{Float64, SparseMatrixCSRtoCSC{T}}(CONE_TYPES)
    end
    function preprocess_constraint(func, set)
        _store_cone_data(dest.cone, set)
        return _preprocess_function(func, set)
    end
    return MOI.copy_to(dest.data, src; preprocess = preprocess_constraint, kws...)
end

function _store_cone_data(::ConeData, ::MOI.Zeros) end
function _store_cone_data(::ConeData, ::MOI.Nonnegatives) end
function _store_cone_data(cone::ConeData, set::MOI.SecondOrderCone)
    push!(cone.qa, set.dimension)
end
function _store_cone_data(cone::ConeData, set::MOI.PositiveSemidefiniteConeTriangle)
    push!(cone.sa, set.side_dimension)
end
function _store_cone_data(::ConeData, ::MOI.ExponentialCone) end
function _store_cone_data(::ConeData, ::MOI.DualExponentialCone) end
function _store_cone_data(cone::ConeData, set::MOI.PowerCone{Float64})
    push!(cone.p, set.exponent)
end
function _store_cone_data(cone::ConeData, set::MOI.DualPowerCone{Float64})
    # SCS' convention: dual cones have a negative exponent.
    push!(cone.p, -set.exponent)
end

# Vectorized length for matrix dimension n
sympackedlen(n) = div(n*(n+1), 2)
# Matrix dimension for vectorized length n
sympackeddim(n) = div(isqrt(1+8n) - 1, 2)
trimap(i::Integer, j::Integer) = i < j ? trimap(j, i) : div((i-1)*i, 2) + j
trimapL(i::Integer, j::Integer, n::Integer) = i < j ? trimapL(j, i, n) : i + div((2n-j) * (j-1), 2)
function _sympackedto(x, n, mapfrom, mapto)
    @assert length(x) == sympackedlen(n)
    y = similar(x)
    for i in 1:n, j in 1:i
        y[mapto(i, j)] = x[mapfrom(i, j)]
    end
    y
end
sympackedLtoU(x, n=sympackeddim(length(x))) = _sympackedto(x, n, (i, j) -> trimapL(i, j, n), trimap)
sympackedUtoL(x, n) = _sympackedto(x, n, trimap, (i, j) -> trimapL(i, j, n))

# Scale coefficients depending on rows index
# rows: List of row indices
# coef: List of corresponding coefficients
# d: dimension of set
# rev: if true, we unscale instead (e.g. divide by √2 instead of multiply for PSD cone)
function _scalecoef(rows::AbstractVector{<: Integer}, coef::Vector{Float64}, d::Integer, rev::Bool)
    scaling = rev ? 1 / √2 : 1 * √2
    output = copy(coef)
    for i in 1:length(output)
        if !MOI.Utilities.is_diagonal_vectorized_index(rows[i])
            output[i] *= scaling
        end
    end
    return output
end

# Unscale the coefficients in `coef` with respective rows in `rows` for a set `s`
scalecoef(rows, coef, s) = _scalecoef(rows, coef, MOI.dimension(s), false)
# Unscale the coefficients in `coef` with respective rows in `rows` for a set of type `S` with dimension `d`
unscalecoef(coef) = _scalecoef(eachindex(coef), coef, sympackeddim(length(coef)), true)

function _scale(i, coef)
    if MOI.Utilities.is_diagonal_vectorized_index(i)
        return coef
    else
        return coef * √2
    end
end

function _preprocess_function(func, set::MOI.PositiveSemidefiniteConeTriangle)
    n = set.side_dimension
    LtoU_map = sympackedLtoU(1:sympackedlen(n), n)
    function map_term(t::MOI.VectorAffineTerm)
        return MOI.VectorAffineTerm(
            LtoU_map[t.output_index],
            MOI.ScalarAffineTerm(
                _scale(t.output_index, t.scalar_term.coefficient),
                t.scalar_term.variable_index
            )
        )
    end
    UtoL_map = sympackedUtoL(1:sympackedlen(n), n)
    function constant(row)
        i = UtoL_map[row]
        return _scale(i, func.constants[i])
    end
    new_func = MOI.VectorAffineFunction{Float64}(
        MOI.VectorAffineTerm{Float64}[map_term(t) for t in func.terms],
        constant.(eachindex(func.constants))
    )
    # The rows have been reordered in `map_term` so we need to re-canonicalize to reorder the rows.
    MOI.Utilities.canonicalize!(new_func)
    return new_func
end
_preprocess_function(func, set) = func

function MOI.optimize!(optimizer::Optimizer)
    data = optimizer.data
    m = data.A.m
    n = data.A.n
    b = data.b
    objective_constant = data.objective_constant
    c = data.c
    if data.sense == MOI.MAX_SENSE
        c = -c
    end

    options = optimizer.options
    if optimizer.silent
        options = copy(options)
        options[:verbose] = 0
    end

    linear_solver, options = sanitize_SCS_options(options)
    P = spzeros(n, n)

    cone = optimizer.cone
    sol = SCS_solve(linear_solver, m, n, _managed_matrix(data.A),
                    ManagedSCSMatrix{Int64}(n, n, P), b, c,
                    data.num_rows[1], data.num_rows[2],
                    Float64[], Float64[],
                    cone.qa, cone.sa, div(data.num_rows[5], 3), div(data.num_rows[6], 3), cone.p,
                    data.primal, data.dual,
                    data.slack; options...)

    data.primal = sol.x
    data.dual = sol.y
    data.slack = sol.s
    ret_val = sol.ret_val
    objective_value = (data.sense == MOI.MAX_SENSE ? -1 : 1) * sol.info.pobj
    dual_objective_value = (data.sense == MOI.MAX_SENSE ? -1 : 1) * sol.info.dobj
    solve_time = (sol.info.setupTime + sol.info.solveTime) / 1000
    optimizer.sol = MOISolution(ret_val, raw_status(sol.info),
                                objective_value, dual_objective_value,
                                objective_constant, solve_time, sol.info.iter)
    return
end

function MOI.get(optimizer::Optimizer, ::MOI.SolveTime)
    return optimizer.sol.solve_time_sec
end
function MOI.get(optimizer::Optimizer, ::MOI.RawStatusString)
    return optimizer.sol.raw_status
end

function MOI.get(optimizer::Optimizer, ::ADMMIterations)
    return optimizer.sol.iterations
end

# Implements getter for result value and statuses
# SCS returns one of the following integers:
# -7 SCS_INFEASIBLE_INACCURATE
# -6 SCS_UNBOUNDED_INACCURATE
# -5 SCS_SIGINT
# -4 SCS_FAILED
# -3 SCS_INDETERMINATE
# -2 SCS_INFEASIBLE  : primal infeasible, dual unbounded
# -1 SCS_UNBOUNDED   : primal unbounded, dual infeasible
#  0 SCS_UNFINISHED  : never returned, used as placeholder
#  1 SCS_SOLVED
#  2 SCS_SOLVED_INACCURATE
function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    s = optimizer.sol.ret_val
    @assert -7 <= s <= 2
    if s == -7
        return MOI.ALMOST_INFEASIBLE
    elseif s == -6
        return MOI.ALMOST_DUAL_INFEASIBLE
    elseif s == 2
        return MOI.ALMOST_OPTIMAL
    elseif s == -5
        return MOI.INTERRUPTED
    elseif s == -4
        return MOI.NUMERICAL_ERROR
    elseif s == -3
        return MOI.SLOW_PROGRESS
    elseif s == -2
        return MOI.INFEASIBLE
    elseif s == -1
        return MOI.DUAL_INFEASIBLE
    elseif s == 1
        return MOI.OPTIMAL
    else
        @assert s == 0
        return MOI.OPTIMIZE_NOT_CALLED
    end
end

function MOI.get(optimizer::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    value = optimizer.sol.objective_value
    if !MOIU.is_ray(MOI.get(optimizer, MOI.PrimalStatus()))
        value += optimizer.sol.objective_constant
    end
    return value
end
function MOI.get(optimizer::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    value = optimizer.sol.dual_objective_value
    if !MOIU.is_ray(MOI.get(optimizer, MOI.DualStatus()))
        value += optimizer.sol.objective_constant
    end
    return value
end

function MOI.get(optimizer::Optimizer, attr::MOI.PrimalStatus)
    if attr.N > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    s = optimizer.sol.ret_val
    if s in (-3, 1, 2)
        MOI.FEASIBLE_POINT
    elseif s in (-6, -1)
        MOI.INFEASIBILITY_CERTIFICATE
    else
        MOI.INFEASIBLE_POINT
    end
end
function MOI.get(optimizer::Optimizer, attr::MOI.VariablePrimal, vi::VI)
    MOI.check_result_index_bounds(optimizer, attr)
    optimizer.data.primal[vi.value]
end
function MOI.get(optimizer::Optimizer, attr::MOI.VariablePrimal, vi::Vector{VI})
    return MOI.get.(optimizer, attr, vi)
end
function post_process_result(x, ::Type{MOI.PositiveSemidefiniteConeTriangle})
    return unscalecoef(sympackedLtoU(x))
end
post_process_result(x, ::Type) = x
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintPrimal, ci::CI{F, S}) where {F, S}
    MOI.check_result_index_bounds(optimizer, attr)
    return post_process_result(optimizer.data.slack[rows(optimizer.data, ci)], S)
end

function MOI.get(optimizer::Optimizer, attr::MOI.DualStatus)
    if attr.N > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    s = optimizer.sol.ret_val
    if s in (-3, 1, 2)
        MOI.FEASIBLE_POINT
    elseif s in (-7, -2)
        MOI.INFEASIBILITY_CERTIFICATE
    else
        MOI.INFEASIBLE_POINT
    end
end
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintDual, ci::CI{F, S}) where {F, S}
    MOI.check_result_index_bounds(optimizer, attr)
    return post_process_result(optimizer.data.dual[rows(optimizer.data, ci)], S)
end

MOI.get(optimizer::Optimizer, ::MOI.ResultCount) = 1
