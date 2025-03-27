module NCTSSoS

using DynamicPolynomials
using DynamicPolynomials:
    AbstractVariable, variables, coefficient, monomial, terms, isconstant, AbstractPolynomial, NonCommutative, Monomial
using SparseArrays, LinearAlgebra
using JuMP
using CliqueTrees
using CliqueTrees: EliminationAlgorithm, SupernodeType
import CliqueTrees.cliquetree
using ChordalGraph
using Graphs
using Base.Iterators: product, flatten

export PolyOpt, StatePolyOpt
export NCStatePolynomial
export TRACE, EIGEN
export SolverConfig
export cs_nctssos, cs_nctssos_higher

# TODO: include other methods in docs
export MF, MMD, NoElimination
export @ncpolyvar, @polyvar

include("statepolynomial.jl")

include("pop.jl")

include("sparse.jl")

include("moment_solver.jl")

include("sos_solver.jl")

include("solver_utils.jl")

include("interface.jl")

end
