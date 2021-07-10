#############################################################################
# SCS.jl
# Wrapper around the SCS solver https://github.com/cvxgrp/scs
#############################################################################
# SCSSolverInterface.jl
# MathProgBase.jl interface for the SCS.jl solver wrapper
#############################################################################

using LinearAlgebra: dot
using MathProgBase.SolverInterface

import MathProgBase.SolverInterface: ConicModel, LinearQuadraticModel,
    getdual, getobjval, getsolution, getsolvetime, getvardual, loadproblem!,
    numconstr, numvar, optimize!, setbvec!, setwarmstart!, status,
    supportedcones

import Base.convert

# TODO: don't add to Base.convert!
function convert(x::Type{Int}, y::UnitRange{Int})
    if length(y) == 1
        return y[1]
    else
        error("convert` has no method matching convert(::Type{Int}, ::UnitRange{Int})")
    end
end
#############################################################################
# Define the MPB Solver and Model objects
export SCSSolver
struct SCSSolver <: AbstractMathProgSolver
    options
end
SCSSolver(;kwargs...) = SCSSolver(kwargs)

mutable struct SCSMathProgModel <: AbstractConicModel
    m::Int                            # Number of constraints
    n::Int                            # Number of variables
    input_numconstr::Int
    input_numvar::Int
    A::SparseMatrixCSC{Float64,Int}   # The A matrix (equalities)
    b::Vector{Float64}                # RHS
    c::Vector{Float64}                # The objective coeffs (always min)
    f::Int                            # number of zero cones
    l::Int                            # number of linear cones { x | x >= 0}
    q::Vector{Int}                    # Array of SOC sizes
    s::Vector{Int}                    # Array of SDP sizes
    ep::Int                           # Number of primal exponential cones
    ed::Int                           # Number of dual exponential cones
    orig_sense::Symbol                # Original objective sense
    # Post-solve
    solve_stat::Symbol
    solve_time::Float64
    obj_val::Float64
    primal_sol::Vector{Float64}
    dual_sol::Vector{Float64}
    slack::Vector{Float64}
    row_map_ind::Vector{Int}
    row_map_type::Vector{Symbol}
    col_map_ind::Vector{Int}          # map from MPB variables to rows
    col_map_type::Vector{Symbol}
    options
end

SCSMathProgModel(;kwargs...) = SCSMathProgModel(0, 0, 0, 0, spzeros(0, 0), Int[], Int[],
                                      0, 0, Int[], Int[], 0, 0,
                                      :Min, :NotSolved, 0.0, 0.0, Float64[], Float64[],
                                      Float64[], Int[], Symbol[],
                                      Int[], Symbol[], kwargs)

#############################################################################
# Begin implementation of the MPB low-level interface
# Implements
# - ConicModel
# - loadproblem!
# - optimize!
# - status
# - numvar
# - numconstr
# http://mathprogbasejl.readthedocs.org/en/latest/solverinterface.html

ConicModel(s::SCSSolver) = SCSMathProgModel(;s.options...)
LinearQuadraticModel(s::SCSSolver) = ConicToLPQPBridge(ConicModel(s))

#=
function setsense!(m::SCSMathProgModel, sns::Symbol)
    if m.orig_sense != sns
        sns == :Min || sns == :Max || error("Unrecognized sense $sns")
        m.orig_sense = sns
        m.c *= -1
    end
    nothing
end
=#

# TODO needs to be updated for newest constants
const status_map = Dict{Int, Symbol}(
    1 => :Optimal,
    -2 => :Infeasible,
    -1 => :Unbounded,
    -3 => :Indeterminate,
    -4 => :Error
)

function optimize!(m::SCSMathProgModel)
    linear_solver, options = sanitize_SCS_options(m.options)
    t = time()
    solution = SCS_solve(linear_solver, m.m, m.n, m.A, m.b, m.c,
                         m.f, m.l, Float64[], Float64[],
                         m.q, m.s, m.ep, m.ed, Float64[],
                         m.primal_sol, m.dual_sol, m.slack; options...)
    m.solve_time = time() - t

    m.solve_stat = get(status_map, solution.ret_val, :UnknownError)
    m.primal_sol = solution.x

    m.dual_sol = solution.y

    # TODO: Get the right slack variables in the right order
    m.slack = solution.s

    m.obj_val = dot(m.c, m.primal_sol) * (m.orig_sense == :Max ? -1 : +1)
end

status(m::SCSMathProgModel) = m.solve_stat
getobjval(m::SCSMathProgModel) = m.obj_val
getsolution(m::SCSMathProgModel) = copy(m.primal_sol)

function invertsdconesize(p)
    return (sqrt(8*p+1) - 1) / 2
end

function isintegertol(n)
    return abs(n - convert(Int, n)) < 1e-16
end

#############################################################################
# Begin implementation of the MPB conic interface
# Implements
# - loadconicproblem!
# - supportedcones
# http://mathprogbasejl.readthedocs.org/en/latest/conic.html

function orderconesforscs(A_in, b_in, c_cones, v_cones)
    # Order the cones as:
    # Free, Zero, NonNeg (NonPos are converted), SOC, SDP, ExpPrimal, ExpDual
    #
    # Returns:
    # - scs_A (A ordered as needed), b
    # - num_free, num_zero, num_linear, soc_sizes, soc_sizes, sqrt_sdp_sizes,
    # - sqrt_sdp_size, num_expprimal, num_expdual

    m, n = size(A_in)
    A_in_t = A_in'
    A_t = spzeros(n,0)
    b = zeros(0)
    row_map_ind = zeros(Int, length(b_in))
    row_map_type = Array{Symbol}(undef, length(b_in))
    col_map_ind = zeros(Int, n)
    col_map_type = Array{Symbol}(undef, n)

    # First, count the total number of variables
    num_vars = 0
    for (cone, idxs) in v_cones
        col_map_type[idxs] .= cone
        num_vars += length(idxs)
    end
    @assert num_vars == n

    num_free = 0
    zeroidx = Int[]
    nonnegidx = Int[]
    nonposidx = Int[]
    socidx = Int[]
    soc_sizes = Int[]
    new_c_cones = Any[]

    for (cone, idxs) in c_cones
        if cone == :Free
            error("Why are you passing in a free constraint?")
        end
        # merge some cones for efficiency
        if cone == :Zero
            append!(zeroidx, idxs)
        elseif cone == :NonNeg
            append!(nonnegidx, idxs)
        elseif cone == :NonPos
            append!(nonposidx, idxs)
        elseif cone == :SOC
            append!(socidx, idxs)
            push!(soc_sizes, length(idxs))
        else
            push!(new_c_cones, (cone,idxs))
        end
    end
    length(zeroidx) > 0 && push!(new_c_cones, (:Zero, zeroidx))
    length(nonnegidx) > 0 && push!(new_c_cones, (:NonNeg, nonnegidx))
    length(nonposidx) > 0 && push!(new_c_cones, (:NonPos, nonposidx))
    length(socidx) > 0 && push!(new_c_cones, (:SOC, socidx))

    for (cone, idxs) in v_cones
        if cone == :Free
            num_free += length(idxs)
        end
    end

    num_zero = 0
    for (cone, idxs) in new_c_cones
        if cone == :Zero
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            A_t = [A_t A_in_t[:,idxs]]
            b = [b; b_in[idxs,:]]
            num_zero += length(idxs)
        end
    end
    for (cone, idxs) in v_cones
        if cone == :Zero
            nidx = length(idxs)
            col_map_ind[idxs] = (length(b)+1):(length(b)+nidx)
            A_t = [A_t sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]
            num_zero += nidx
        end
    end

    num_lin = 0
    for (cone, idxs) in new_c_cones
        if cone == :NonNeg
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            A_t = [A_t A_in_t[:,idxs]]
            b = [b; b_in[idxs,:]]
            num_lin += length(idxs)
        elseif cone == :NonPos
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            A_t = [A_t -A_in_t[:,idxs]]
            b = [b; -b_in[idxs,:]]
            num_lin += length(idxs)
        end
    end
    for (cone, idxs) in v_cones
        nidx = length(idxs)
        col_map_ind[idxs] = (length(b)+1):(length(b)+nidx)
        if cone == :NonNeg
            A_t = [A_t -sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]
            num_lin += nidx
        elseif cone == :NonPos
            A_t = [A_t sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]
            num_lin += nidx
        end
    end

    for (cone, idxs) in new_c_cones
        if cone == :SOC
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            A_t = [A_t A_in_t[:,idxs]]
            b = [b; b_in[idxs,:]]
        end
    end
    for (cone, idxs) in v_cones
        if cone == :SOC
            nidx = length(idxs)
            col_map_ind[idxs] = (length(b)+1):(length(b)+nidx)
            A_t = [A_t -sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]
            push!(soc_sizes, nidx)
        end
    end

    sqrt_sdp_sizes = Int[]
    for (cone, idxs) in new_c_cones
        if cone == :SDP
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            A_t = [A_t A_in_t[:,idxs]]
            b = [b; b_in[idxs,:]]
            # n must be a square integer
            n = length(idxs)
            isintegertol(invertsdconesize(n)) || error("number of SDP variables must be n*(n+1)/2")
            sqrt_n = convert(Int, invertsdconesize(n));
            push!(sqrt_sdp_sizes, sqrt_n)
        end
    end
    for (cone, idxs) in v_cones
        if cone == :SDP
            nidx = length(idxs)
            col_map_ind[idxs] = (length(b)+1):(length(b)+nidx)
            A_t = [A_t -sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]
             # n must be a square integer
            isintegertol(invertsdconesize(nidx)) || error("number of SDP variables must be n*(n+1)/2")
            sqrt_n = convert(Int, invertsdconesize(nidx));
            push!(sqrt_sdp_sizes, sqrt_n)
        end
    end

    num_expprimal = 0
    for (cone, idxs) in new_c_cones
        if cone == :ExpPrimal
            length(idxs) % 3 == 0 ||
                error("Number of ExpPrimal variables must be a multiple of 3")
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            A_t = [A_t A_in_t[:,idxs]]
            b = [b; b_in[idxs,:]]

            num_expprimal += div(length(idxs), 3)
        end
    end
    for (cone, idxs) in v_cones
        if cone == :ExpPrimal
            length(idxs) % 3 == 0 ||
                error("Number of ExpPrimal variables must be a multiple of 3")
            nidx = length(idxs)
            col_map_ind[idxs] = (length(b)+1):(length(b)+nidx)
            A_t = [A_t -sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]

            num_expprimal += div(length(idxs), 3)
        end
    end

    num_expdual = 0
    for (cone, idxs) in new_c_cones
        if cone == :ExpDual
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            length(idxs) % 3 == 0 ||
                error("Number of ExpDual variables must be a multiple of 3")
            A_t = [A_t A_in_t[:,idxs]]
            b = [b; b_in[idxs,:]]

            num_expdual += div(length(idxs), 3)
        end
    end
    for (cone, idxs) in v_cones
        if cone == :ExpDual
            length(idxs) % 3 == 0 ||
                error("Number of ExpDual variables must be a multiple of 3")
            nidx = length(idxs)
            col_map_ind[idxs] = (length(b)+1):(length(b)+nidx)
            A_t = [A_t -sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]

            num_expdual += div(length(idxs), 3)
        end
    end

    return A_t', b, num_free, num_zero, num_lin, soc_sizes,
           sqrt_sdp_sizes, num_expprimal, num_expdual, col_map_ind, col_map_type, row_map_ind, row_map_type
end


loadproblem!(model::SCSMathProgModel, c, A, b, constr_cones, var_cones) =
    loadproblem!(model, c, sparse(A), b, constr_cones, var_cones)


function loadproblem!(model::SCSMathProgModel, c, A::SparseMatrixCSC, b, constr_cones, var_cones)
    # TODO: We should support SOCRotated
    bad_cones = [:SOCRotated]
    for cone_vars in constr_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) not supported")
    end
    for cone_vars in var_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) not supported")
    end

    # Convert idxs to an array
    c_cones = [(cone, [idxs...]) for (cone, idxs) in constr_cones]
    v_cones = [(cone, [idxs...]) for (cone, idxs) in var_cones]

    scs_A, scs_b, num_free, f, l, q, s, ep, ed, col_map_ind, col_map_type, row_map_ind, row_map_type =
        orderconesforscs(A, b, c_cones, v_cones)

    m, n = size(scs_A)

    model.n             = n
    model.m             = m # + rows_G
    model.A             = scs_A
    model.b             = scs_b[:]
    model.c             = c[:]
    model.q             = q
    model.s             = s
    model.ep            = ep
    model.ed            = ed
    model.orig_sense    = :Min
    model.f             = f
    model.l             = l
    model.col_map_ind   = col_map_ind
    model.col_map_type  = col_map_type
    model.row_map_ind   = row_map_ind
    model.row_map_type  = row_map_type
    model.input_numconstr = size(A,1)
    model.input_numvar    = size(A,2)

    return model
end

numvar(model::SCSMathProgModel) = model.input_numvar
numconstr(model::SCSMathProgModel) = model.input_numconstr

supportedcones(s::SCSSolver) = [:Free, :Zero, :NonNeg, :NonPos, :SOC, :SDP, :ExpPrimal, :ExpDual]

function getdual(m::SCSMathProgModel)
    dual = m.dual_sol[m.row_map_ind]
    # flip sign for NonPos since it's treated as NonNeg by SCS
    for i in 1:length(m.row_map_type)
        if m.row_map_type[i] == :NonPos
            dual[i] = -dual[i]
        end
    end
    return dual
end

function getvardual(m::SCSMathProgModel)
    dual = zeros(length(m.col_map_ind))
    for i in 1:length(m.col_map_type)
        if m.col_map_type[i] == :Free
            continue # dual is zero
        elseif m.col_map_type[i] == :NonPos
            # flip sign for NonPos since it's treated as NonNeg by SCS
            dual[i] = -m.dual_sol[m.col_map_ind[i]]
        else
            dual[i] = m.dual_sol[m.col_map_ind[i]]
        end
    end
    return dual
end

function addoption!(m::SCSMathProgModel, option::Symbol, value)
    nt = NamedTuple{(option,), Tuple{typeof(value)}}((value,))
    m.options = pairs(merge(m.options.data, nt))
    return m
end

# warmstart
# kwargs can be `primal_sol`, `dual_sol`, and `slack`
function setwarmstart!(m::SCSMathProgModel, primal_sol; kwargs...)
    addoption!(m, :warm_start, true)
    m.primal_sol = primal_sol
    for (k,v) in kwargs
        setfield!(m, k, v)
    end

    # check sizes to prevent segfaults
    nconstr, nvar = size(m.A)
    length(m.primal_sol) == nvar || (m.primal_sol = zeros(nvar))
    length(m.dual_sol) == nconstr || (m.dual_sol = zeros(nconstr))
    length(m.slack) == nconstr || (m.slack = zeros(nconstr))
    m
end

function setbvec!(m::SCSMathProgModel, b::Vector{Float64})
    m.b[m.row_map_ind] = b
end

getsolvetime(m::SCSMathProgModel) = m.solve_time
