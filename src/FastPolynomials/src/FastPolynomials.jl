module FastPolynomials

export variables, ς

export @ncpolyvar, Variable
include("variables.jl")

export Monomial
include("monomials.jl")

export Polynomial
include("polynomial.jl")

include("compare.jl")

include("arithmetic.jl")

include("utils.jl")

include("state_word.jl")

include("statepolynomial.jl")

end
