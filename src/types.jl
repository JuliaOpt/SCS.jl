using SparseArrays

export SCSMatrix, SCSData, SCSSettings, SCSSolution, SCSInfo, SCSCone

VecOrMatOrSparse = Union{VecOrMat, SparseMatrixCSC{Float64,Int}}

SCSInt = Union{Int32, Int64}
abstract type LinearSolver end
struct DirectSolver <: LinearSolver end
struct IndirectSolver <: LinearSolver end
struct IndirectGpuSolver <: LinearSolver end

scsint_t(::Type{<:LinearSolver}) = Int
scsint_t(::Type{IndirectGpuSolver}) = Int32

clib(::Type{DirectSolver}) = direct
clib(::Type{IndirectSolver}) = indirect
clib(::Type{IndirectGpuSolver}) = indirectgpu


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

    function ManagedSCSMatrix{T}(m::Integer, n::Integer, values::Vector{Cdouble},
                                 rowval::Vector{T}, colptr::Vector{T}) where T
        # scsmatref holds the reference to SCSMatrix created out of data in ManagedSCSMatrix.
        # this way the reference to SCSMatrix always points to valid data
        # as long as the ManagedSCSMatrix is not GC collected.
        # One MUST
        #   `Base.unsafe_convert(Ref{SCSMatrix}, scsmatref)`
        # when assigning to a field of type `Ptr{SCSMatrix}`, or specify
        #   `Ref{SCSMatrix}`
        # in the type tuple when ccalling with `scsmatref`

        scsmat = SCSMatrix{T}(
            pointer(values), pointer(rowval), pointer(colptr), m, n)

        return new{T}(values, rowval, colptr, Base.cconvert(Ref{SCSMatrix{T}}, scsmat))
    end
end

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

function ManagedSCSMatrix{T}(m::Integer, n::Integer, A::AbstractMatrix) where T
    return ManagedSCSMatrix{T}(m, n, sparse(A))
end

struct SCSSettings{T<:SCSInt}
    normalize::T # boolean, heuristic data rescaling
    scale::Cdouble # if normalized, rescales by this factor
    rho_x::Cdouble # x equality constraint scaling
    max_iters::T # maximum iterations to take
    eps::Cdouble # convergence tolerance
    alpha::Cdouble # relaxation parameter
    cg_rate::Cdouble # for indirect, tolerance goes down like (1/iter)^cg_rate
    verbose::T # boolean, write out progress
    warm_start::T # boolean, warm start (put initial guess in Sol struct)
    acceleration_lookback::T # acceleration memory parameter
    adaptive_scaling::T # whether to adaptively update the scale param
    write_data_filename::Cstring
    log_csv_filename::Cstring

    SCSSettings{T}() where T = new{T}()
    function SCSSettings{T}(normalize, scale, rho_x, max_iters, eps, alpha,
                            cg_rate, verbose, warm_start,
                            acceleration_lookback, adaptive_scaling,
                            write_data_filename, log_csv_filename) where T
        return new{T}(normalize, scale, rho_x, max_iters, eps, alpha, cg_rate,
                      verbose, warm_start, acceleration_lookback, adaptive_scaling,
                      write_data_filename, log_csv_filename)
    end
end

function _SCS_user_settings(default_settings::SCSSettings{T};
        normalize=default_settings.normalize,
        scale=default_settings.scale,
        rho_x=default_settings.rho_x,
        max_iters=default_settings.max_iters,
        eps=default_settings.eps,
        alpha=default_settings.alpha,
        cg_rate=default_settings.cg_rate,
        verbose=default_settings.verbose,
        warm_start=default_settings.warm_start,
        acceleration_lookback=default_settings.acceleration_lookback,
        adaptive_scaling=default_settings.adaptive_scaling,
        write_data_filename=C_NULL,
        log_csv_filename=C_NULL
        ) where T
    write_data_filename_p = write_data_filename isa String ? pointer(write_data_filename) : C_NULL
    log_csv_filename_p = log_csv_filename isa String ? pointer(log_csv_filename) : C_NULL
    return SCSSettings{T}(normalize, scale, rho_x, max_iters, eps, alpha,
                         cg_rate, verbose,warm_start, acceleration_lookback,
                         adaptive_scaling, write_data_filename_p,
                         log_csv_filename_p)
end

# Warning: Strings provided through options must outlive the solve.
# SCSSettings keeps only a pointer to the strings.
function SCSSettings(linear_solver::Type{<:LinearSolver}; options...)

    T = scsint_t(linear_solver)

    default_settings = Base.cconvert(Ref{SCSSettings{T}}, SCSSettings{T}())
    dummy_data = SCSData{T}(0,0,C_NULL, C_NULL, C_NULL, C_NULL,
        Base.unsafe_convert(Ref{SCSSettings{T}}, default_settings))

    Base.GC.@preserve default_settings begin
        SCS_set_default_settings(linear_solver, dummy_data)
    end

    return _SCS_user_settings(default_settings[]; options...)
end

struct SCSData{T<:SCSInt}
    # A has m rows, n cols
    m::T
    n::T
    A::Ptr{SCSMatrix{T}}
    P::Ptr{SCSMatrix{T}}
    # b is of size m, c is of size n
    b::Ptr{Cdouble}
    c::Ptr{Cdouble}
    stgs::Ptr{SCSSettings{T}}
end

function SCSData(m::Integer, n::Integer,
    A::ManagedSCSMatrix{T},
    P::ManagedSCSMatrix{T},
    b::Vector{Float64},
    c::Vector{Float64},
    stgs::Ref{SCSSettings{T}}) where T
    return SCSData{T}(m, n,
        Base.unsafe_convert(Ref{SCSMatrix{T}}, A.scsmatref), # Ptr{SCSMatrix{T}}
        Base.unsafe_convert(Ref{SCSMatrix{T}}, P.scsmatref), # Ptr{SCSMatrix{T}}
        pointer(b),
        pointer(c),
        Base.unsafe_convert(Ref{SCSSettings{T}}, stgs) # Ptr{SCSSettings{T}}
    )
end

struct SCSSolution
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
    s::Ptr{Cdouble}
end

struct SCSInfo{T<:SCSInt}
    iter::T
    status::NTuple{64, Cchar} # char status[64]
    statusVal::T
    pobj::Cdouble
    dobj::Cdouble
    resPri::Cdouble
    resDual::Cdouble
    resInfeas::Cdouble
    resUnbdd::Cdouble
    xt_p_x::Cdouble
    relGap::Cdouble
    setupTime::Cdouble
    solveTime::Cdouble
end

SCSInfo{T}() where T = SCSInfo{T}(0, ntuple(_ -> zero(Cchar), 64), 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

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
# box cones { (t, x) | t * l <= x <= t * u, t >= 0 },
# second-order cones { (t,x) | ||x||_2 <= t },
# semi-definite cones { X | X psd }, and
# exponential cones {(x,y,z) | y e^(x/y) <= z, y>0 }.
# dual exponential cones {(u,v,w) | −u e^(v/u) <= e w, u<0}
# power cones {(x,y,z) | x^a * y^(1-a) >= |z|, x>=0, y>=0}
# dual power cones {(u,v,w) | (u/a)^a * (v/(1-a))^(1-a) >= |w|, u>=0, v>=0}
struct SCSCone{T<:SCSInt}
    f::T # number of linear equality constraints
    l::T # length of LP cone
    bu::Ptr{Cdouble} # upper box values
    bl::Ptr{Cdouble} # lower box values
    bsize::T # length of box constraint, including scale t
    q::Ptr{T} # array of second-order cone constraints
    qsize::T # length of SOC array
    s::Ptr{T} # array of SD constraints
    ssize::T # length of SD array
    ep::T # number of primal exponential cone triples
    ed::T # number of dual exponential cone triples
    p::Ptr{Cdouble} # array of power cone params, must be \in [-1, 1], negative values are interpreted as specifying the dual cone
    psize::T # length of power cone array
end

# Returns an SCSCone. The q, s, and p arrays are *not* GC tracked in the
# struct. Use this only when you know that q, s, and p will outlive the struct.
function SCSCone{T}(f::Integer, l::Integer, bu::Vector{Cdouble},
                    bl::Vector{Cdouble}, q::Vector{T}, s::Vector{T},
                    ep::Integer, ed::Integer, p::Vector{Cdouble}) where T
    @assert length(bl) == length(bu)
    return SCSCone{T}(f, l, pointer(bu), pointer(bl),
                      isempty(bu) ? 0 : length(bu) + 1,
                      pointer(q), length(q), pointer(s), length(s), ep, ed,
                      pointer(p), length(p))
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
        linear_solver = IndirectSolver # the default linear_solver
    end

    SCS_options = append!([:linear_solver], fieldnames(SCSSettings))
    unrecognized = setdiff(keys(options), SCS_options)
    if length(unrecognized) > 0
        plur = length(unrecognized) > 1 ? "s" : ""
        throw(ArgumentError("Unrecognized option$plur passed to SCS: $(join(unrecognized, ", "));\nRecognized options are: $(join(SCS_options, ", ", " and "))."))
    end
    return linear_solver, options
end
