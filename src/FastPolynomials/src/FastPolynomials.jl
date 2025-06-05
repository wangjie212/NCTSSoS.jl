module FastPolynomials

export @ncpolyvar, Variable
include("variables.jl")

export Monomial
include("monomials.jl")

export Polynomial
include("polynomial.jl")

include("compare.jl")


include("arithmetic.jl")

include("operators.jl")

end
