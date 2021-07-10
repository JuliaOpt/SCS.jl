@testset "Test the MathProgBase wrapper with linprog" begin
    include("mpb_linear.jl")
end

import MathProgBase
for T in solvers
    @testset "MathProgBase $T" begin
        include(joinpath(dirname(dirname(pathof(MathProgBase))), "test", "conicinterface.jl"))
        coniclineartest(SCS.SCSSolver(linear_solver=T, eps_abs=1e-6, eps_rel=1e-6, verbose=0),
                        duals=true, tol=1e-5)
        conicSOCtest(SCS.SCSSolver(linear_solver=T, eps_abs=1e-6, eps_rel=1e-6, verbose=0),
                     duals=true, tol=1e-5)
        conicEXPtest(SCS.SCSSolver(linear_solver=T, eps_abs=1e-6, eps_rel=1e-7, verbose=0),
                     duals=true, tol=1e-5)
        conicSDPtest(SCS.SCSSolver(linear_solver=T, eps_abs=1e-6, eps_rel=1e-7, verbose=0),
                     duals=true, tol=1e-5)
    end
end
