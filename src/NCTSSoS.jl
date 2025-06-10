module NCTSSoS

using SparseArrays, LinearAlgebra, JuMP
using CliqueTrees, ChordalGraph, Graphs
using CliqueTrees: EliminationAlgorithm, SupernodeType
import CliqueTrees.cliquetree

include("FastPolynomials/src/FastPolynomials.jl")
using .FastPolynomials
using .FastPolynomials: Variable, Monomial, Polynomial, StateWord, NCStateWord, StatePolynomial, NCStatePolynomial
using .FastPolynomials: sorted_union, monomials, _comm, sorted_unique, _projective, _unipotent
using .FastPolynomials: monomials, maxdegree, get_basis, symmetric_canonicalize, neat_dot
using .FastPolynomials: monomials, coefficients, terms, get_state_basis
export @ncpolyvar


export PolyOpt, TRACE, EIGEN, StatePolyOpt
export SolverConfig
export NoElimination, MF, MMD
export cs_nctssos

include("pop.jl")

include("solver_utils.jl")

include("sparse.jl")

include("moment_solver.jl")

include("sos_solver.jl")

include("interface.jl")

include("state_sparse.jl")

include("state_moment_solver.jl")

include("state_sos_solver.jl")


end
