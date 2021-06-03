function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    # TODO: Base.precompile(Tuple{Core.kwftype(typeof(copy_to)),NamedTuple{(:copy_names,), Tuple{Bool}},typeof(MathOptInterface.copy_to),Optimizer,MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}})   # time: 0.383372
    Base.precompile(Tuple{typeof(MathOptInterface.optimize!),Optimizer})   # time: 0.17350855
    Base.precompile(Tuple{typeof(_allocate_constraints),GeometricConicForm{Float64, SparseMatrixCSRtoCSC{Int64}, NTuple{8, DataType}},MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},MathOptInterface.Utilities.IndexMap,Int64,Type{MathOptInterface.Zeros}})   # time: 0.11354891
    Base.precompile(Tuple{typeof(_allocate_constraints),GeometricConicForm{Float64, SparseMatrixCSRtoCSC{Int64}, NTuple{8, DataType}},MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},MathOptInterface.Utilities.IndexMap,Int64,Type{MathOptInterface.PositiveSemidefiniteConeTriangle}})   # time: 0.10451437
    Base.precompile(Tuple{typeof(_allocate_constraints),GeometricConicForm{Float64, SparseMatrixCSRtoCSC{Int64}, NTuple{8, DataType}},MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},MathOptInterface.Utilities.IndexMap,Int64,Type{MathOptInterface.Nonnegatives}})   # time: 0.06301721
    Base.precompile(Tuple{Core.kwftype(typeof(SCS_solve)),NamedTuple{(:verbose,), Tuple{Int64}},typeof(SCS_solve),Type{IndirectSolver},Int64,Int64,ManagedSCSMatrix{Int64},Vector{Float64},Vector{Float64},Int64,Int64,Vector{Int64},Vector{Int64},Int64,Int64,Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64}})   # time: 0.051085986
    Base.precompile(Tuple{typeof(MathOptInterface.get),Optimizer,MathOptInterface.ConstraintPrimal,MathOptInterface.ConstraintIndex{_A, _B} where {_A, _B}})   # time: 0.048825916
    Base.precompile(Tuple{typeof(MathOptInterface.empty!),Optimizer})   # time: 0.021038031
    Base.precompile(Tuple{typeof(_allocate_constraints),GeometricConicForm{Float64, SparseMatrixCSRtoCSC{Int64}, NTuple{8, DataType}},MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}},MathOptInterface.Utilities.IndexMap,Int64,Type{MathOptInterface.SecondOrderCone}})   # time: 0.012096747
    Base.precompile(Tuple{var"##copy_to#38",Function,Bool,typeof(MathOptInterface.copy_to),GeometricConicForm{Float64, SparseMatrixCSRtoCSC{Int64}, NTuple{8, DataType}},MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}})   # time: 0.011085093
    Base.precompile(Tuple{typeof(MathOptInterface.supports_constraint),Optimizer,Type{MathOptInterface.VectorAffineFunction{Float64}},Type{MathOptInterface.Nonpositives}})   # time: 0.006092956
    Base.precompile(Tuple{typeof(MathOptInterface.supports_constraint),Optimizer,Type{MathOptInterface.VectorAffineFunction{Float64}},Type{Union{}}})   # time: 0.005302456
    Base.precompile(Tuple{typeof(MathOptInterface.supports_constraint),Optimizer,Type{MathOptInterface.VectorAffineFunction{Float64}},Type{MathOptInterface.PositiveSemidefiniteConeSquare}})   # time: 0.005300446
    Base.precompile(Tuple{typeof(MathOptInterface.supports_constraint),Optimizer,Type{MathOptInterface.VectorAffineFunction{Float64}},Type{MathOptInterface.Reals}})   # time: 0.005227482
    Base.precompile(Tuple{typeof(MathOptInterface.get),Optimizer,MathOptInterface.PrimalStatus})   # time: 0.004503549
    Base.precompile(Tuple{typeof(MathOptInterface.get),Optimizer,MathOptInterface.ObjectiveValue})   # time: 0.004384848
    Base.precompile(Tuple{typeof(MathOptInterface.set),Optimizer,MathOptInterface.Silent,Bool})   # time: 0.003145233
    Base.precompile(Tuple{typeof(MathOptInterface.set),GeometricConicForm{Float64, SparseMatrixCSRtoCSC{Int64}, NTuple{8, DataType}},MathOptInterface.ObjectiveFunction{MathOptInterface.ScalarAffineFunction{Float64}},MathOptInterface.ScalarAffineFunction{Float64}})   # time: 0.00264804
    Base.precompile(Tuple{typeof(MathOptInterface.supports_constraint),Optimizer,Type{MathOptInterface.VectorAffineFunction{Float64}},Type{MathOptInterface.Zeros}})   # time: 0.0023151
    Base.precompile(Tuple{typeof(MathOptInterface.supports_constraint),Optimizer,Type{MathOptInterface.VectorAffineFunction{Float64}},Type{MathOptInterface.Nonnegatives}})   # time: 0.002270974
    Base.precompile(Tuple{typeof(MathOptInterface.supports_constraint),Optimizer,Type{MathOptInterface.VectorAffineFunction{Float64}},Type{MathOptInterface.PositiveSemidefiniteConeTriangle}})   # time: 0.002066653
    Base.precompile(Tuple{typeof(MathOptInterface.supports_constraint),Optimizer,Type{MathOptInterface.VectorAffineFunction{Float64}},Type{MathOptInterface.SecondOrderCone}})   # time: 0.001853893
    Base.precompile(Tuple{typeof(_managed_matrix),SparseMatrixCSRtoCSC{Int64}})   # time: 0.001720187
    Base.precompile(Tuple{typeof(MathOptInterface.get),Optimizer,MathOptInterface.VariablePrimal,MathOptInterface.VariableIndex})   # time: 0.001075106
end
