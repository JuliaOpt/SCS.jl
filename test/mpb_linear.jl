#############################################################################
# SCS.jl
# Wrapper around the SCS solver https://github.com/cvxgrp/scs
# See http://github.com/karanveerm/SCS.jl
#############################################################################
# test/mpb_linear.jl
# Test the MathProgBase.jl interface for the SCS.jl solver wrapper
#############################################################################

using MathProgBase
using LinearAlgebra
using SparseArrays

solver = SCSSolver(eps_abs=1e-5, eps_rel=1e-5)
objtol = 1e-4
primaltol = 1e-4

# min -x
# s.t. 2x + y <= 1.5
# x,y >= 0
# solution is (0.75,0) with objval -0.75
sol = linprog([-1,0],[2 1],'<',1.5,solver)
@test sol.status == :Optimal
@test isapprox(sol.objval, -0.75, atol=objtol)
@test isapprox(norm(sol.sol - [0.75,0.0]), 0, atol=primaltol)

sol = linprog([-1,0],sparse([2 1]),'<',1.5,solver)
@test sol.status == :Optimal
@test isapprox(sol.objval, -0.75, atol=objtol)
@test isapprox(norm(sol.sol - [0.75,0.0]), 0, atol=primaltol)

# test infeasible problem:
# min x
# s.t. 2x+y <= -1
# x,y >= 0
sol = linprog([1.0,0.0],[2.0 1.0],'<',-1.0,solver)
@test sol.status == :Infeasible

# test unbounded problem:
# min -x-y
# s.t. -x+2y <= 0
# x,y >= 0
sol = linprog([-1.,-1.],[-1. 2.],'<',[0.0],solver)
@test sol.status == :Unbounded

# More complicated through loadproblem!
A = [1.0 1.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 1.0;
     0.0 0.0 1.0 1.0 1.0]
collb = [0.0, 0.0, 0.0, 0.0, 0.0]
colub = [Inf, Inf, Inf, Inf, Inf]
obj   = [3.0, 4.0, 4.0, 9.0, 5.0]
sense = :Max
rowlb = [-Inf, -Inf, -Inf]
rowub = [ 5.0,  3.0,  9.0]
m = MathProgBase.LinearQuadraticModel(solver)
MathProgBase.loadproblem!(m, A, collb, colub, obj, rowlb, rowub, sense)
MathProgBase.optimize!(m)
@test isapprox(MathProgBase.getobjval(m), 99.0, atol=1e-3)
x = MathProgBase.getsolution(m)
@test isapprox(x[1], 2.0, atol=1e-4)
@test isapprox(x[2], 3.0, atol=1e-4)
@test isapprox(x[3], 0.0, atol=1e-4)
@test isapprox(x[4], 9.0, atol=1e-4)
@test isapprox(x[5], 0.0, atol=1e-4)
