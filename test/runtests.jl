using Test
using SCS

solvers = SCS.available_solvers

include("test_problems.jl")

function feasible_basic_problems(solver)
    A = reshape([1.0],(1,1))
    P = spzeros(1, 1)
    solution = SCS_solve(solver, 1, 1, A, P, [1.0], [1.0], 1, 0, Float64[], Float64[], Int[], Int[], 0, 0, Float64[])
    @test solution.ret_val == 1
    feasible_basic_conic(solver)
    feasible_exponential_conic(solver)
    feasible_sdp_conic(solver)
    feasible_pow_conic(solver)
end

for s in solvers
    feasible_basic_problems(s)
end

#include("options.jl")

include("MOI_wrapper.jl")
#include("MPB_wrapper.jl")
