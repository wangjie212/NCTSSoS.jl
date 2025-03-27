struct StateWord{V,M}
    state_monos::Vector{Monomial{V,M}}
end

Base.show(io::IO, sw::StateWord) = print(io, join(map(x -> "⟨$(x)⟩", sw.state_monos), " "))

struct StatePolynomial{V,M,T}
    coeffs::Vector{T}
    state_word::Vector{StateWord{V,M}}
end

Base.show(io::IO, sp::StatePolynomial) = print(io, join(map(x -> "$(x[1]) * $(x[2])", zip(sp.coeffs, sp.state_word)), " + "))

# T: type of coefficient
# V: whether the variabels is NonCommutative{CreationOrder} or Commutative{CreationOrder}
# M: ordering of the monomials Graded{LexOrder} or else
# FIXME: inheriting from AbstractPolynomial gives StackOverflowError when trying to print in REPL
struct NCStatePolynomial{V,M,T}
    state_poly::Vector{StatePolynomial{V,M,T}}
    words::Vector{Monomial{V,M}}

    # TODO: get a more user friendly way of declaring StatePolynomial
    function NCStatePolynomial(coeffs::Vector{Vector{T}}, formal_words::Vector{Vector{Vector{Monomial{V,M}}}}, words::Vector{Monomial{V,M}}) where {V,M,T}
        @assert length(coeffs) == length(formal_words) == length(words) "Coefficients, formal words, and words must have the same length"
        new{V,M,T}(coeffs, formal_words, words)
    end
end

function NCStatePolynomial(formal_part::Vector{Vector{Polynomial{V,M,T}}}, words::Vector{Monomial{V,M}}) where {V,M,T}
    coeffs, formal_words = mapfoldl((args...) -> push!.(args...), formal_part, init=(Vector{T}[], Vector{Vector{Monomial{V,M}}}[])) do formal_poly_vec
        terms_tuple_vec = reshape(collect(product(map(x -> terms(x), formal_poly_vec)...)), :)
        coeff_vec, monovec = mapfoldl((args...) -> push!.(args...), terms_tuple_vec; init=(T[], Vector{Monomial{V,M}}[])) do terms_tuple
            reduce(*, coefficient.(terms_tuple)), [monomial(t) for t in terms_tuple]
        end
        coeff_vec, monovec
    end
    @show coeffs formal_words words
    @show length(coeffs) length(formal_words) length(words)
    return StatePolynomial(coeffs, formal_words, words)
end

Base.show(io::IO, sp::NCStatePolynomial) = print(io, join(map(zip(sp.coeffs, sp.formal_words, sp.words)) do (coeff, formal_word, word)
        formal_word_str = "(" * join(map(zip(coeff, formal_word)) do (cc, ffs)
                "$(cc) ⋅ " * join(map(x -> "⟨$(string(x))⟩", ffs), " ⋅ ")
            end, " + ")
        formal_word_str * ") ⋅ $(string(word))"
    end, " + "))


DynamicPolynomials.variables(sp::NCStatePolynomial) = sorted_union(variables.(sp.words)..., map(x -> union(map(y -> union(variables.(y)...), x)...), sp.formal_words)...)
