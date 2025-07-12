module FastPolynomials

using MLStyle: @match

export variables, ς, tr

export @ncpolyvar, Variable
include("variables.jl")

export Monomial, monomial
include("monomials.jl")

export Polynomial, degree, monomials, coefficients
include("polynomial.jl")

include("compare.jl")

include("arithmetic.jl")

include("utils.jl")

include("phased_monomial.jl")

include("state_word.jl")

export ncstatepoly
include("statepolynomial.jl")

export SimplifyAlgorithm, simplify, canonicalize
include("simplify.jl")
end
