# T: type of coefficient
# V: whether the variabels is NonCommutative{CreationOrder} or Commutative{CreationOrder}
# M: ordering of the monomials Graded{LexOrder} or else
# FIXME: inheriting from AbstractPolynomial gives StackOverflowError when trying to print in REPL
struct StatePolynomial{V,M,T} 
    coeffs::Vector{T}
    formal_words::Vector{Vector{Monomial{V,M}}} # expectation value of operators w.r.t general state
    words::Vector{Monomial{V,M}}
end

# TODO: get a more user friendly way of declaring StatePolynomial
function StatePolynomial(coeffs::Vector{T}, formal_words::Vector{Vector{Monomial{V,M}}}, words::Vector{Monomial{V,M}}) where {V,M,T}
    @assert length(coeffs) == length(formal_words) == length(words) "Coefficients, formal words, and words must have the same length"
    StatePolynomial{V,M,T}(coeffs, formal_words, words)
end

Base.show(io::IO, sp::StatePolynomial) = print(io, join(map(zip(sp.coeffs,sp.formal_words,sp.words)) do (coeff, formal_word, word)
    formal_word_str = join(map(x->"⟨$(x)⟩",string.(formal_word)), " ⋅ ")
   "$(coeff) ⋅ $(formal_word_str) ⋅ $(word)" 
end, " + "))