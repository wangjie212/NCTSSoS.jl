"""
   PolynomialLike

Type alias for objects that can be either a Variable or a Monomial or a Polynomial.
Used in arithmetic operations to allow flexible input types.
"""
const PolynomialLike = Union{Variable,Monomial,Polynomial}

function Base.:(+)(a::PolynomialLike, b::PolynomialLike)
    ap = Polynomial(a)
    bp = Polynomial(b)
    return Polynomial([ap.coeffs; bp.coeffs], [ap.monos; bp.monos])
end

Base.:(+)(a::PolynomialLike,b::Number) = a + Polynomial(b)
Base.:(+)(a::Number,b::PolynomialLike) = Polynomial(a) + b


function Base.:(-)(a::PolynomialLike, b::PolynomialLike)
    ap = Polynomial(a)
    bp = Polynomial(b)
    return Polynomial([ap.coeffs; -bp.coeffs], [ap.monos; bp.monos])
end

Base.:(-)(a::PolynomialLike,b::Number) = a - Polynomial(b)
Base.:(-)(a::Number,b::PolynomialLike) = Polynomial(a) - b

function Base.:(-)(a::Polynomial)
    return Polynomial(-a.coeffs, a.monos)
end

"""
    Base.:(*)(a::Polynomial, b::Polynomial)

Multiplies two polynomials using distributive property.

# Arguments
- `a::Polynomial`: First polynomial
- `b::Polynomial`: Second polynomial

# Returns
- `Polynomial`: Product polynomial with all pairwise coefficient and monomial products
"""
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


"""
    Base.:(*)(a::Polynomial{T}, b::PolynomialLike) where {T}

Multiplies a polynomial by a variable/monomial.

# Arguments
- `a::Polynomial{T}`: Polynomial
- `b::PolynomialLike`: Variable or monomial

# Returns
- `Polynomial{T}`: Product with coefficient type T
"""
Base.:(*)(a::PolynomialLike, b::PolynomialLike) = *(promote(a, b)...)

Base.:(*)(a::Variable, b::Variable) = Monomial([a, b], [1, 1])

function Base.:(*)(x::Monomial, y::Monomial)
    i = findlast(z -> z > 0, x.z)
    isnothing(i) && return y
    j = findfirst(z -> z > 0, y.z)
    isnothing(j) && return x

    if x.vars[i] == y.vars[j]
        w = [x.vars[1:i]; y.vars[(j + 1):end]]
        z = [x.z[1:(i - 1)]; x.z[i] + y.z[j]; y.z[(j + 1):end]]
    else
        w = [x.vars[1:i]; y.vars[j:end]]
        z = [x.z[1:i]; y.z[j:end]]
    end
    return Monomial(w, z)
end

Base.:(*)(a::Number, b::Polynomial) = Polynomial(a .* b.coeffs, b.monos)
Base.:(*)(a::Polynomial, b::Number) = Polynomial(b .* a.coeffs, a.monos)
Base.:(*)(a::Number, b::PolynomialLike) = a * Polynomial(b)
Base.:(*)(a::PolynomialLike, b::Number) = Polynomial(a) * b

Base.:(/)(a::PolynomialLike, b::Number) = Polynomial(a) * inv(b)
