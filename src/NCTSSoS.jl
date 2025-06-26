module NCTSSoS

using SparseArrays, LinearAlgebra, JuMP
using CliqueTrees, ChordalGraph, Graphs
using CliqueTrees: EliminationAlgorithm, SupernodeType
import CliqueTrees.cliquetree

include("FastPolynomials/src/FastPolynomials.jl")
using .FastPolynomials
using .FastPolynomials: AbstractPolynomial
using .FastPolynomials:
    Variable,
    Monomial,
    Polynomial,
    StateWord,
    NCStateWord,
    StatePolynomial,
    NCStatePolynomial
using .FastPolynomials:
    sorted_union, monomials, _comm, sorted_unique, _projective, _unipotent
using .FastPolynomials: monomials, maxdegree, get_basis, symmetric_canonicalize, neat_dot
using .FastPolynomials: monomials, coefficients, terms, get_state_basis
using .FastPolynomials: expval
export @ncpolyvar


# export PolyOpt, TRACE, EIGEN
# export SolverConfig
# export NoElimination, MF, MMD, AsIsElimination, MaximalElimination
# export cs_nctssos, cs_nctssos_higher

# include("pop.jl")

# include("elimination.jl")

# include("solver_utils.jl")

# include("sparse.jl")

# include("moment_solver.jl")

# include("sos_solver.jl")

# include("interface.jl")
end
