if get(ENV, "BUILDKITE", "false") == "true"
    # This file requires a GPU in order to run. It gets tested as part of the
    # JuliaGPU CI on Buildkite. Contact @odow for more details.
    include("test_gpu.jl")
    exit(0)
end

using Test
using SCS

solvers = SCS.available_solvers

include("test_problems.jl")
include("MOI_wrapper.jl")

@testset "SCS" begin
    @testset "Basic feasible problems: $s" for s in solvers
        feasible_basic_problems(s)
    end

    include("options.jl")
    include("MPB_wrapper.jl")

    @testset "MOI wrapper: $s" for s in solvers
        moi_tests(s)
    end
end
