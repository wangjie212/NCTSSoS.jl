module NCTSSoS

using JuMP

include("FastPolynomials/src/FastPolynomials.jl")
using .FastPolynomials
export @ncpolyvar


export PolyOpt, TRACE, EIGEN

include("pop.jl")

include("solver_utils.jl")


end
