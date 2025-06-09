module NCTSSoS

using JuMP
using CliqueTrees
using CliqueTrees: EliminationAlgorithm, SupernodeType
import CliqueTrees.cliquetree

include("FastPolynomials/src/FastPolynomials.jl")
using .FastPolynomials
export @ncpolyvar


export PolyOpt, TRACE, EIGEN

include("pop.jl")

include("solver_utils.jl")

include("moment_solver.jl")


end
