module NCTSSoS

include("FastPolynomials/src/FastPolynomials.jl")
using .FastPolynomials
export @ncpolyvar

<<<<<<< HEAD
<<<<<<< HEAD
export PolyOpt, TRACE, EIGEN
<<<<<<< HEAD
=======
=======
=======

include("pop.jl")

include("solver_utils.jl")

>>>>>>> 42f8261 (fix conflict)
# using DynamicPolynomials
# using DynamicPolynomials:
#     AbstractVariable, variables, coefficient, monomial, terms, isconstant, AbstractPolynomial, NonCommutative, Monomial, degree
# using SparseArrays, LinearAlgebra
# using JuMP
# using CliqueTrees
# using CliqueTrees: EliminationAlgorithm, SupernodeType
# import CliqueTrees.cliquetree
# using ChordalGraph
# using Graphs
# using Base.Iterators: product, flatten
# using DynamicPolynomials.MP: compare

<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 220ab62 (rebase and get docs)
export PolyOpt, StatePolyOpt
export NCStatePolynomial, ς, NCStateWord
export TRACE, EIGEN
export SolverConfig
export cs_nctssos, cs_nctssos_higher
=======
# export PolyOpt, StatePolyOpt
# export NCStatePolynomial, ς
# export TRACE, EIGEN
# export SolverConfig
# export cs_nctssos, cs_nctssos_higher
>>>>>>> 31cac3a (up)

=======
>>>>>>> 42f8261 (fix conflict)
# # TODO: include other methods in docs
# export MF, MMD, NoElimination, MinimalChordal
# export @ncpolyvar, @polyvar

# # Export polynomial types
# export NCStateWord, NCStateTerm
# export StatePolynomial, NCStatePolynomial
# export PolyOpt, StatePolyOpt

# # Export sparsity types
# export CorrelativeSparsity, StateCorrelativeSparsity
# export TermSparsity, StateTermSparsity

# # Export solver types
# export MomentProblem, StateMomentProblem
# export SOSProblem, PolyOptResult

# # Export operations
# export symmetric_canonicalize, cyclic_canonicalize
# export neat_dot, get_basis, support
# export reducer, expval

<<<<<<< HEAD
include("statepolynomial.jl")
>>>>>>> f3d3314 (update State Polynomial Example)
=======
# include("statepolynomial.jl")
>>>>>>> 220ab62 (rebase and get docs)

# include("pop.jl")

<<<<<<< HEAD
include("solver_utils.jl")

# using DynamicPolynomials
# using DynamicPolynomials:
#     AbstractVariable, variables, coefficient, monomial, terms, isconstant, AbstractPolynomial, NonCommutative, Monomial, degree
# using SparseArrays, LinearAlgebra
# using JuMP
# using CliqueTrees
# using CliqueTrees: EliminationAlgorithm, SupernodeType
# import CliqueTrees.cliquetree
# using ChordalGraph
# using Graphs
# using Base.Iterators: product, flatten
# using DynamicPolynomials.MP: compare

# export StatePolyOpt
# export NCStatePolynomial, ς
# export SolverConfig
# export cs_nctssos, cs_nctssos_higher

# # TODO: include other methods in docs
# export MF, MMD, NoElimination, MinimalChordal
# export @ncpolyvar, @polyvar

# # Export polynomial types
# export NCStateWord, NCStateTerm
# export StatePolynomial, NCStatePolynomial
# export PolyOpt, StatePolyOpt

# # Export sparsity types
# export CorrelativeSparsity, StateCorrelativeSparsity
# export TermSparsity, StateTermSparsity

# # Export solver types
# export MomentProblem, StateMomentProblem
# export SOSProblem, PolyOptResult

# # Export operations
# export symmetric_canonicalize, cyclic_canonicalize
# export neat_dot, get_basis, support
# export reducer, expval

# include("statepolynomial.jl")


# include("sparse.jl")

# include("state_sparse.jl")

# include("moment_solver.jl")

# include("state_moment_solver.jl")

# include("sos_solver.jl")

# include("state_sos_solver.jl")

<<<<<<< HEAD

=======
# include("sparse.jl")

# include("state_sparse.jl")

# include("moment_solver.jl")

# include("state_moment_solver.jl")

# include("sos_solver.jl")

# include("state_sos_solver.jl")

# include("solver_utils.jl")
=======
>>>>>>> 42f8261 (fix conflict)

>>>>>>> 220ab62 (rebase and get docs)
# include("interface.jl")

end
