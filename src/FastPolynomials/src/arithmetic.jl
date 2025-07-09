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

function _concat_var_expos(
    a::Vector{Variable}, a_z::Vector{Int}, b::Vector{Variable}, b_z::Vector{Int}
)
    i = findlast(z -> z > 0, a_z)
    isnothing(i) && return (b, b_z)
    j = findfirst(z -> z > 0, b_z)
    isnothing(j) && return (a, a_z)

    lb = length(b)
    if a[i] == b[j]
        total_len = i + lb - j
        w = Vector{Variable}(undef, total_len)
        z = Vector{Int}(undef, total_len)
        @inbounds for k in 1:(i - 1)
            w[k] = a[k]
            z[k] = a_z[k]
        end

        w[i] = a[i]
        z[i] = a_z[i] + b_z[j]

        @inbounds for k in (j + 1):lb
            w[i + k - j] = b[k]
            z[i + k - j] = b_z[k]
        end
        # w = [view(a, 1:i); view(b, (j + 1):lb)]
        # z = [view(a_z, 1:(i - 1)); a_z[i] + b_z[j]; view(b_z, (j + 1):lb)]
    else
        total_len = i + lb - j + 1
        w = Vector{Variable}(undef, total_len)
        z = Vector{Int}(undef, total_len)
        @inbounds for k in 1:i
            w[k] = a[k]
            z[k] = a_z[k]
        end
        @inbounds for k in j:lb
            w[i + k - j + 1] = b[k]
            z[i + k - j + 1] = b_z[k]
        end
        # w = [view(a, 1:i); view(b, j:length(b))]
        # z = [view(a_z, 1:i); view(b_z, j:length(b_z))]
    end
    return (w, z)
end

function Base.:(*)(x::Monomial, y::Monomial)
    w, z = _concat_var_expos(x.vars, x.z, y.vars, y.z)
    return Monomial(w, z)
end

function Base.:(*)(x::AdjointMonomial, y::Monomial)
    w, z = _concat_var_expos(reverse(x.parent.vars), reverse(x.parent.z), y.vars, y.z)
    return Monomial(w, z)
end

Base.:(*)(a::Number, b::Polynomial) = Polynomial(a .* b.coeffs, b.monos)
Base.:(*)(a::Polynomial, b::Number) = Polynomial(b .* a.coeffs, a.monos)
Base.:(*)(a::Number, b::PolynomialLike) = a * Polynomial(b)
Base.:(*)(a::PolynomialLike, b::Number) = Polynomial(a) * b

Base.:(/)(a::PolynomialLike, b::Number) = Polynomial(a) * inv(b)
