module FastPolynomials

export variables, ς

export @ncpolyvar, Variable
include("variables.jl")

export Monomial, monomial
include("monomials.jl")

export Polynomial, degree
include("polynomial.jl")

include("compare.jl")

include("arithmetic.jl")

include("utils.jl")

include("state_word.jl")

export ncstatepoly
include("statepolynomial.jl")

end
