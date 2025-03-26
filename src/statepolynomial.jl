# TODO: update to conform with DynamicPolynomial v0.5.0 and up
# T: type of coefficient
# V: whether the variabels is commutative
# M: ordering of the monomials Lexico or Reverse (this occurs in v0.5.0 of DynamicPolynomials)
struct StatePolynomial{V,T} <: AbstractPolynomial{T}
    coeffs::Vector{T}
    formal_words::Vector{MonomialVector{V}} # expectation value of operators w.r.t general state
    words::MonomialVector{V}
end

# TODO: get a more user friendly way of declaring StatePolynomial
function StatePolynomial{V,T}(coeffs::Vector{T}, formal_words::Vector{Vector{Monomial{V}}}, words::Vector{Monomial{V}}) where {V,T}
    @assert length(coeffs) == length(formal_words) == length(words) "Coefficients, formal words, and words must have the same length"
    StatePolynomial{V,T}(coeffs, formal_words, words)
end
