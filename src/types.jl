using SparseArrays

export SCSMatrix, SCSData, SCSSettings, SCSSolution, SCSInfo, SCSCone

VecOrMatOrSparse = Union{VecOrMat, SparseMatrixCSC{Float64,Int}}

SCSInt = Union{Int32, Int64}
abstract type LinearSolver end
struct DirectSolver <: LinearSolver end
struct IndirectSolver <: LinearSolver end
struct GpuIndirectSolver <: LinearSolver end

scsint_t(::Type{<:LinearSolver}) = Int
scsint_t(::Type{GpuIndirectSolver}) = Int32

clib(::Type{DirectSolver}) = direct
clib(::Type{IndirectSolver}) = indirect
clib(::Type{GpuIndirectSolver}) = gpuindirect


struct SCSMatrix{T<:SCSInt}
    values::Ptr{Cdouble}
    rowval::Ptr{T}
    colptr::Ptr{T}
    m::T
    n::T
end

# Version where Julia manages the memory for the vectors.
struct ManagedSCSMatrix{T<:SCSInt}
    values::Vector{Cdouble}
    rowval::Vector{T}
    colptr::Vector{T}
    scsmatref::Base.RefValue{SCSMatrix{T}}

    function ManagedSCSMatrix{T}(
        m::Integer,
        n::Integer,
        values::Vector{Cdouble},
        rowval::Vector{T},
        colptr::Vector{T},
    ) where {T}
        # scsmatref holds the reference to SCSMatrix created out of data in ManagedSCSMatrix.
        # this way the reference to SCSMatrix always points to valid data
        # as long as the ManagedSCSMatrix is not GC collected.
        # One MUST
        #   `Base.unsafe_convert(Ref{SCSMatrix}, scsmatref)`
        # when assigning to a field of type `Ptr{SCSMatrix}`, or specify
        #   `Ref{SCSMatrix}`
        # in the type tuple when ccalling with `scsmatref`

        scsmat = SCSMatrix{T}(pointer(values), pointer(rowval), pointer(colptr), m, n)

        return new{T}(values, rowval, colptr, Base.cconvert(Ref{SCSMatrix{T}}, scsmat))
    end
end
csize(m::ManagedSCSMatrix) = (M = m.scsmatref[]; (M.m, M.n))

function ManagedSCSMatrix{T}(m::Integer, n::Integer, A::SparseMatrixCSC) where T
    # we convert everything to make sure that no conversion will
    # accidentally happen when ccalling `new`, so that `scsmatref`
    # holds pointers to actual data stored in the fields.

    # A reference to A.nzval is kept if no conversion is needed.
    values = convert(Vector{Cdouble}, A.nzval)
    rowval = convert(Vector{T}, A.rowval .- 1)
    colptr = convert(Vector{T}, A.colptr .- 1)

    return ManagedSCSMatrix{T}(m, n, values, rowval, colptr)
end

function ManagedSCSMatrix{T}(m::Integer, n::Integer, A::AbstractMatrix) where {T}
    return ManagedSCSMatrix{T}(m, n, sparse(A))
end

"""
    SCSSettings
SCS struct with parameters controlling scs.

NOTE: The `write_data_filename` and `log_csv_filename` are stored as **pointers** in the struct,
so be careful to ensure that a Julia references to these arrays exists as long as `SCSSettings` will be used.
"""
struct SCSSettings{T<:SCSInt}
    normalize::T # boolean, heuristic data rescaling
    scale::Cdouble # if normalized, set factor for rescales
    rho_x::Cdouble # x equality constraint scaling
    max_iters::T # maximum iterations to take
    eps_abs::Cdouble # absolute convergence tolerance
    eps_rel::Cdouble # relative convergence tolerance
    eps_infeas::Cdouble # infeasible convergence tolerance
    alpha::Cdouble # relaxation parameter
    time_limit_secs::Cdouble # time limit in secs
    verbose::T # boolean, write out progress
    warm_start::T # boolean, warm start (put initial guess in SCSSolution struct)
    acceleration_lookback::T # acceleration memory parameter
    acceleration_interval::T # interval to apply acceleration
    adaptive_scaling::T # whether to adaptively update the scale param
    write_data_filename::Cstring # if set scs will dump raw prob data
    log_csv_filename::Cstring # if set scs will log solve

    SCSSettings{T}() where {T} = new{T}()
    SCSSettings{T}(
        normalize,
        scale,
        rho_x,
        max_iters,
        eps_abs,
        eps_rel,
        eps_infeas,
        alpha,
        time_limit_secs,
        verbose,
        warm_start,
        acceleration_lookback,
        acceleration_interval,
        adaptive_scaling,
        write_data_filename::Union{String, Cstring},
        log_csv_filename::Union{String, Cstring},
    ) where {T} = new{T}(
        normalize,
        scale,
        rho_x,
        max_iters,
        eps_abs,
        eps_rel,
        eps_infeas,
        alpha,
        time_limit_secs,
        verbose,
        warm_start,
        acceleration_lookback,
        acceleration_interval,
        adaptive_scaling,
        Base.unsafe_convert(Cstring, write_data_filename),
        Base.unsafe_convert(Cstring, log_csv_filename),
    )
end

function _SCS_user_settings(
    default_settings::SCSSettings{T};
    normalize = default_settings.normalize,
    scale = default_settings.scale,
    rho_x = default_settings.rho_x,
    max_iters = default_settings.max_iters,
    eps_abs = default_settings.eps_abs,
    eps_rel = default_settings.eps_rel,
    eps_infeas = default_settings.eps_infeas,
    alpha = default_settings.alpha,
    time_limit_secs = default_settings.time_limit_secs,
    verbose = default_settings.verbose,
    warm_start = default_settings.warm_start,
    acceleration_lookback = default_settings.acceleration_lookback,
    acceleration_interval = default_settings.acceleration_interval,
    adaptive_scaling = default_settings.adaptive_scaling,
    write_data_filename = default_settings.write_data_filename,
    log_csv_filename = default_settings.log_csv_filename,
) where {T}
    return SCSSettings{T}(
        normalize,
        scale,
        rho_x,
        max_iters,
        eps_abs,
        eps_rel,
        eps_infeas,
        alpha,
        time_limit_secs,
        verbose,
        warm_start,
        acceleration_lookback,
        acceleration_interval,
        adaptive_scaling,
        write_data_filename,
        log_csv_filename,
    )
end

function SCSSettings(linear_solver::Type{<:LinearSolver}; options...)
    T = scsint_t(linear_solver)
    default_settings = Base.cconvert(Ref{SCSSettings{T}}, SCSSettings{T}())
    scsm_ptr = convert(Ptr{SCSMatrix{T}}, C_NULL)
    ptr_d = convert(Ptr{Cdouble}, C_NULL)
    data = SCSData{T}(0, 0, scsm_ptr, scsm_ptr, ptr_d, ptr_d,
        Base.unsafe_convert(Ref{SCSSettings{T}}, default_settings))

    SCS_set_default_settings(linear_solver, data)

    return _SCS_user_settings(default_settings[]; options...)
end

struct SCSData{T<:SCSInt}
    m::T
    n::T
    A::Ptr{SCSMatrix{T}} # size: (m, n)
    P::Ptr{SCSMatrix{T}} # size: (n, n)
    b::Ptr{Cdouble} # size: (m,1)
    c::Ptr{Cdouble} # size: (n,1)
    stgs::Ptr{SCSSettings{T}}
end

"""
    SCSData
SCS struct with problem definition.

NOTE: The `A`, `P`, `b`, `c` and `stgs` are stored as **pointers** in the struct, so be careful to ensure
that a Julia references to these arrays exists as long as `SCSData` will be used.
"""
function SCSData(
    m::Integer,
    n::Integer,
    A::ManagedSCSMatrix{T},
    P::ManagedSCSMatrix{T},
    b::Vector{Float64},
    c::Vector{Float64},
    stgs::Ref{SCSSettings{T}},
) where {T}
    @assert csize(A) == (size(b, 1), size(c, 1)) == (m, n)
    @assert csize(P) == (n, n)
    return SCSData{T}(
        m,
        n,
        Base.unsafe_convert(Ref{SCSMatrix{T}}, A.scsmatref), # Ptr{SCSMatrix{T}}
        Base.unsafe_convert(Ref{SCSMatrix{T}}, P.scsmatref), # Ptr{SCSMatrix{T}}
        pointer(b),
        pointer(c),
        Base.unsafe_convert(Ref{SCSSettings{T}}, stgs), # Ptr{SCSSettings{T}}
    )
end

struct SCSSolution
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
    s::Ptr{Cdouble}
end

struct SCSInfo{T<:SCSInt}
    iter::T # number of iterations taken
    status::NTuple{64,Cchar} # status string, e.g. 'solved'
    statusVal::T # status as scs_int, defined in scs/inlcude/glbopts.h
    scale_updates::T # number of updates to scale
    pobj::Cdouble # primal objective
    dobj::Cdouble # dual objective
    resPri::Cdouble # primal equality residual
    resDual::Cdouble # dual equality residual
    resInfeas::Cdouble # infeasibility cert residual
    resUnbdd_a::Cdouble # unbounded cert residual
    resUnbdd_p::Cdouble # unbounded cert residual
    relGap::Cdouble # relative duality gap
    setupTime::Cdouble # time taken for setup phase (milliseconds)
    solveTime::Cdouble # time taken for solve phase (milliseconds)
    scale::Cdouble # (final) scale parameter
end

SCSInfo{T}() where {T} = SCSInfo{T}(
    0,
    ntuple(_ -> zero(Cchar), 64),
    0,
    0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
)

function raw_status(info::SCSInfo)
    s = collect(info.status)
    len = findfirst(iszero, s) - 1
    # There is no method String(::Vector{Cchar}) so we convert to `UInt8`.
    return String(UInt8[s[i] for i in 1:len])
end

# SCS solves a problem of the form
# minimize        c' * x
# subject to      A * x + s = b
#                 s in K
# where K is a product cone of
# zero cones,
# linear cones { x | x >= 0 },
# second-order cones { (t,x) | ||x||_2 <= t },
# semi-definite cones { X | X psd }, and
# exponential cones {(x,y,z) | y e^(x/y) <= z, y>0 }.
# dual exponential cones {(u,v,w) | −u e^(v/u) <= e w, u<0}
# power cones {(x,y,z) | x^a * y^(1-a) >= |z|, x>=0, y>=0}
# dual power cones {(u,v,w) | (u/a)^a * (v/(1-a))^(1-a) >= |w|, u>=0, v>=0}
struct SCSCone{T<:SCSInt}
    f::T # number of linear equality constraints
    l::T # length of LP cone
    bu::Ptr{Cdouble} # upper box values, length = bsize
    bl::Ptr{Cdouble} # lower box values, length = bsize
    bsize::T # length of box cone constraint, including scale t
    q::Ptr{T} # array of second-order cone constraints
    qsize::T # length of SOC array
    s::Ptr{T} # array of SD constraints
    ssize::T # length of SD array
    ep::T # number of primal exponential cone triples
    ed::T # number of dual exponential cone triples
    p::Ptr{Cdouble} # array of power cone params, must be ∈ [-1, 1], negative values are interpreted as specifying the dual cone
    psize::T # number of (primal and dual) power cone triples
end

"""
    SCSCone{T}
SCS struct with cone data.

NOTE: The `bu`, `bl`, `q`, `s`, and `p` are stored as **pointers** in the struct,
so be careful to ensure that a Julia references to these arrays exists as long as `SCSCone` will be used.
"""
function SCSCone{T}(
    f::Integer,
    l::Integer,
    bu::Vector{Cdouble},
    bl::Vector{Cdouble},
    q::Vector{T},
    s::Vector{T},
    ep::Integer,
    ed::Integer,
    p::Vector{Cdouble},
) where {T}
    @assert length(bu) == length(bl)
    return SCSCone{T}(
        f,
        l,
        pointer(bu),
        pointer(bl),
        length(bu),
        pointer(q),
        length(q),
        pointer(s),
        length(s),
        ep,
        ed,
        pointer(p),
        length(p),
    )
end


mutable struct Solution{T<:SCSInt}
    x::Array{Float64, 1}
    y::Array{Float64, 1}
    s::Array{Float64, 1}
    info::SCSInfo{T}
    ret_val::T
end

function sanitize_SCS_options(options)
    options = Dict(options)
    if haskey(options, :linear_solver)
        linear_solver = options[:linear_solver]
        if !(linear_solver in available_solvers)
            throw(ArgumentError("Unrecognized linear_solver passed to SCS: $linear_solver;\nRecognized options are: $(join(available_solvers, ", ", " and "))."))
        end
        delete!(options, :linear_solver)
    else
        linear_solver = DirectSolver # the default linear_solver
    end

    SCS_options = append!([:linear_solver], fieldnames(SCSSettings))
    unrecognized = setdiff(keys(options), SCS_options)
    if length(unrecognized) > 0
        plur = length(unrecognized) > 1 ? "s" : ""
        throw(ArgumentError("Unrecognized option$plur passed to SCS: $(join(unrecognized, ", "));\nRecognized options are: $(join(SCS_options, ", ", " and "))."))
    end
    return linear_solver, options
end
