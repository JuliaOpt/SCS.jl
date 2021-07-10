# SCS solves a problem of the form
# minimize        c' * x
# subject to      A * x + s = b
#                 s in K
# where K is a product cone of
# zero cone { x | x = 0 },
# positive orthant { x | x >= 0 },
# second-order cones (SOC) { (t,x) | ||x||_2 <= t },
# semi-definite cones (SDC) { X | X psd }, and
# exponential cones {(x,y,z) | y e^(x/y) <= z, y>0 },
# power cone { (x,y,z) | x^a*y^(1-a) >= |z|, x>=0, y>=0 },
# dual power cone { (u,v,w) | (u/a)^a * (v/(1-a))^(1-a) >= |w|, u>=0, v>=0}
#
#
# Description of input argments:
# A is the matrix with m rows and n cols
# b is a vector of length m
# c is a vector of length n
#
# The rows of A correspond to cones in K and need to be specified in the order above.
#
# f (num primal zero / dual free cones, i.e. primal equality constraints)
# l (num linear cones)
# bu (box cone constraints upper)
# bl (box cone constraints lower)
# q (array of SOCs sizes)
# s (array of SDCs sizes)
# ep (num primal exponential cones)
# ed (num dual exponential cones)
# p (array of power cone params, must be in [-1, 1], negative values specify the dual cone)
#
# Returns a Solution object.
function SCS_solve(linear_solver::Type{<:LinearSolver},
        m::Integer,
        n::Integer,
        A::ManagedSCSMatrix{T},
        P::ManagedSCSMatrix{T},
        b::Vector{Float64},
        c::Vector{Float64},
        f::Integer,
        l::Integer,
        bu::Vector{Float64},
        bl::Vector{Float64},
        q::Vector{<:Integer},
        s::Vector{<:Integer},
        ep::Integer,
        ed::Integer,
        p::Vector{Float64},
        primal_sol::Vector{Float64}=zeros(n),
        dual_sol::Vector{Float64}=zeros(m),
        slack::Vector{Float64}=zeros(m);
        options...) where T

    n > 0 || throw(ArgumentError("The number of variables in SCSModel must be greater than 0"))
    m > 0 || throw(ArgumentError("The number of constraints in SCSModel must be greater than 0"))

    ws = (:warm_start=>true) in options
    warmstart_sizes_are_correct = length(primal_sol) == n && length(dual_sol) == length(slack) == m

    if warmstart_sizes_are_correct
        if !ws
            fill!(primal_sol, 0.0)
            fill!(dual_sol, 0.0)
            fill!(slack, 0.0)
        end
    else
        ws && throw(ArgumentError("The provided warmstart doesn't match problem sizes"))
        primal_sol = zeros(n)
        dual_sol = zeros(m)
        slack = zeros(m)
    end

    @assert T == scsint_t(linear_solver)
    q_T = convert(Vector{T}, q)
    s_T = convert(Vector{T}, s)

    solution = SCSSolution(pointer(primal_sol), pointer(dual_sol), pointer(slack))

    settings = Base.cconvert(Ref{SCSSettings{T}}, SCSSettings(linear_solver; options...))
    # settings potentially hold pointers to strings from options that need to be protected from GC

    data = SCSData(m, n, A, P, b, c, settings)
    # data holds pointers to objects which need to be protected from GC:
    # A, P, b, c and settings

    cone = SCSCone{T}(f, l, bu, bl, q_T, s_T, ep, ed, p)
    # cone holds pointers to objects which need to be protected from GC:
    # bu, bl, q_T, s_T, p

    info_ref = Base.cconvert(Ref{SCSInfo{T}}, SCSInfo{T}())

    Base.GC.@preserve A P b c settings bu bl q_T s_T p options begin
        p_work = SCS_init(linear_solver, data, cone, info_ref)
        status = SCS_solve(linear_solver, p_work, data, cone, solution, info_ref)
        SCS_finish(linear_solver, p_work)
    end

    return Solution(primal_sol, dual_sol, slack, info_ref[], status)
end
function SCS_solve(linear_solver::Type{<:LinearSolver},
    m::Integer,
    n::Integer,
    A::VecOrMatOrSparse,
    P::VecOrMatOrSparse,
    b::Vector{Float64},
    c::Vector{Float64},
    f::Integer,
    l::Integer,
    bu::Vector{Float64},
    bl::Vector{Float64},
    q::Vector{<:Integer},
    s::Vector{<:Integer},
    ep::Integer,
    ed::Integer,
    p::Vector{Float64},
    primal_sol::Vector{Float64}=zeros(n),
    dual_sol::Vector{Float64}=zeros(m),
    slack::Vector{Float64}=zeros(m);
    options...
)

    T = scsint_t(linear_solver)
    return SCS_solve(linear_solver, m, n,
        ManagedSCSMatrix{T}(m, n, A),
        ManagedSCSMatrix{T}(n, n, P),
        b, c, f, l, bu, bl, q, s, ep, ed, p, primal_sol, dual_sol, slack; options...)
end

function SCS_solve(linear_solver::Type{<:LinearSolver},
    m::Integer,
    n::Integer,
    A::ManagedSCSMatrix{T},
    args...;
    options...) where T
    P = ManagedSCSMatrix{T}(n, n, spzeros(n,n))
    return SCS_solve(linear_solver, m, n, A, P, args...; options...)
end

function SCS_solve(linear_solver::Type{<:LinearSolver},
    m::Integer,
    n::Integer,
    A::AbstractMatrix,
    args...;
    options...) where T
    mA = ManagedSCSMatrix{scsint_t(linear_solver)}(size(A)..., A)
    return SCS_solve(linear_solver, m, n, mA, args...; options...)
end

# Wrappers for the direct C API.
# Do not call these wrapper methods directly unless you understand the
# use of GC.@preserve in the SCS_solve helper above.

for linear_solver in available_solvers
    lib = clib(linear_solver)
    T = scsint_t(linear_solver)
    @eval begin
        function SCS_set_default_settings(::Type{$linear_solver}, data::SCSData{$T})
            ccall((:scs_set_default_settings, $lib), Cvoid, (Ref{SCSData{$T}},),
            data)
        end

        # data and cone are const in :scs_init
        function SCS_init(::Type{$linear_solver}, data::SCSData{$T}, cone::SCSCone{$T}, info_ref::Ref{SCSInfo{$T}})

            p_work = ccall((:scs_init, $lib), Ptr{Cvoid},
                (Ref{SCSData{$T}}, Ref{SCSCone{$T}}, Ref{SCSInfo{$T}}),
                data, cone, info_ref)

            return p_work
        end

        # data and cone are const in :scs_solve
        # solution struct contains only `Ptr`s so passing by value
        function SCS_solve(::Type{$linear_solver}, p_work::Ptr{Nothing}, data::SCSData{$T}, cone::SCSCone{$T}, solution::SCSSolution, info_ref::Ref{SCSInfo{$T}})

            status = ccall((:scs_solve, $lib), $T,
                (Ptr{Cvoid}, Ref{SCSData{$T}}, Ref{SCSCone{$T}}, Ref{SCSSolution}, Ref{SCSInfo{$T}}),
                p_work, data, cone, solution, info_ref)

            return status
        end

        function SCS_finish(::Type{$linear_solver}, p_work::Ptr{Nothing})
            ccall((:scs_finish, $lib), Cvoid,
                (Ptr{Cvoid}, ),
                p_work)
        end
    end
end

# This one is safe to call
function SCS_version()
    return unsafe_string(ccall((:scs_version, SCS.direct), Cstring, ()))
end
