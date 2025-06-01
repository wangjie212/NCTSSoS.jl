module FastPolynomials

export @ncpolyvar
include("variables.jl")

export Monomial
include("monomials.jl")

export Polynomial
include("polynomial.jl")

include("compare.jl")


include("arithmetic.jl")

end
