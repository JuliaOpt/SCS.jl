using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

# UniversalFallback is needed for starting values
const CACHE = MOIU.UniversalFallback(MOIU.Model{Float64}())

import SCS

for T in solvers
    optimizer = SCS.Optimizer(linear_solver=T, eps=1e-7)
    MOI.set(optimizer, MOI.Silent(), true)

    @testset "SolverName" begin
        @test MOI.get(optimizer, MOI.SolverName()) == "SCS"
    end

    MOI.empty!(CACHE)
    cached = MOIU.CachingOptimizer(CACHE, optimizer)
    bridged = MOIB.full_bridge_optimizer(cached, Float64)
    config = MOIT.TestConfig(atol=1e-5)

    @testset "Unit with $T" begin
        MOIT.unittest(bridged, config, [
            # FIXME `NumberOfThreads` not supported.
            "number_threads",
            # `TimeLimitSec` not supported.
            "time_limit_sec",
            # ArgumentError: The number of constraints in SCSModel must be greater than 0
            "solve_unbounded_model",
            # Integer and ZeroOne sets are not supported
            "solve_integer_edge_cases", "solve_objbound_edge_cases",
            "solve_zero_one_with_bounds_1",
            "solve_zero_one_with_bounds_2",
            "solve_zero_one_with_bounds_3"])
    end

    @testset "Continuous linear problems with $T" begin
        MOIT.contlineartest(bridged, config)
    end

    @testset "ADMMIterations attribute with $T" begin
        MOIT.linear1test(bridged, config)
        @test MOI.get(bridged, SCS.ADMMIterations()) > 0
    end

    @testset "Continuous quadratic problems with $T" begin
        MOIT.qcptest(bridged, config)
    end

    @testset "Continuous conic problems with $T" begin
        MOIT.contconictest(bridged, config, ["rootdets", "logdets"])
    end
end

@testset "MOI.RawParameter" begin
    model = SCS.Optimizer()
    # TODO(odow): remove symbol cases when deprecation is removed.
    MOI.set(model, MOI.RawParameter(:eps), 1.0)
    @test MOI.get(model, MOI.RawParameter(:eps)) == 1.0
    @test MOI.get(model, MOI.RawParameter("eps")) == 1.0
    MOI.set(model, MOI.RawParameter("eps"), 2.0)
    @test MOI.get(model, MOI.RawParameter(:eps)) == 2.0
    @test MOI.get(model, MOI.RawParameter("eps")) == 2.0
end
