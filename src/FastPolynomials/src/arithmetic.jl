const PolynomialLike = Union{Variable,Monomial,Polynomial}

function Base.:(+)(a::PolynomialLike, b::PolynomialLike)
    ap = Polynomial(a)
    bp = Polynomial(b)
    return Polynomial([ap.coeffs; bp.coeffs], [ap.monos; bp.monos])
end

Base.:(+)(a::PolynomialLike, b::Number) = a + Polynomial(b)
Base.:(+)(a::Number, b::PolynomialLike) = Polynomial(a) + b

function Base.:(-)(a::PolynomialLike, b::PolynomialLike)
    ap = Polynomial(a)
    bp = Polynomial(b)
    return Polynomial([ap.coeffs; -bp.coeffs], [ap.monos; bp.monos])
end

Base.:(-)(a::PolynomialLike, b::Number) = a - Polynomial(b)
Base.:(-)(a::Number, b::PolynomialLike) = Polynomial(a) - b

function Base.:(-)(a::Polynomial)
    return Polynomial(-a.coeffs, a.monos)
end

function Base.:(*)(a::Polynomial, b::Polynomial)
    return Polynomial(
        vec([ca * cb for (ca, cb) in Iterators.product(a.coeffs, b.coeffs)]),
        vec([ma * mb for (ma, mb) in Iterators.product(a.monos, b.monos)]),
    )
end

Base.promote_rule(::Type{Monomial}, ::Type{Variable}) = Monomial
Base.promote_rule(::Type{Polynomial{T}}, ::Type{Monomial}) where {T} = Polynomial{T}
Base.promote_rule(::Type{Variable}, ::Type{Polynomial{T}}) where {T} = Polynomial{T}
Base.promote_rule(::Type{Polynomial{T}}, ::Type{Variable}) where {T} = Polynomial{T}

Base.:(*)(a::PolynomialLike, b::PolynomialLike) = *(promote(a, b)...)
Base.:(*)(a::Variable, b::Variable) = monomial([a, b], [1, 1])

function Base.:(*)(x::Monomial, y::Monomial)
    i = findlast(z -> z > 0, x.z)
    isnothing(i) && return y
    j = findfirst(z -> z > 0, y.z)
    isnothing(j) && return x

    if x.vars[i] == y.vars[j]
        w = [view(x.vars, 1:i); view(y.vars, (j + 1):length(y.vars))]
        z = [view(x.z, 1:(i - 1)); x.z[i] + y.z[j]; view(y.z, (j + 1):length(y.z))]
    else
        w = [view(x.vars, 1:i); view(y.vars, j:length(y.vars))]
        z = [view(x.z, 1:i); view(y.z, j:length(y.z))]
    end
    return monomial(w, z)
end

Base.:(*)(a::Number, b::Polynomial) = Polynomial(a .* b.coeffs, b.monos)
Base.:(*)(a::Polynomial, b::Number) = Polynomial(b .* a.coeffs, a.monos)
Base.:(*)(a::Number, b::PolynomialLike) = a * Polynomial(b)
Base.:(*)(a::PolynomialLike, b::Number) = Polynomial(a) * b

Base.:(/)(a::PolynomialLike, b::Number) = Polynomial(a) * inv(b)
