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
    isempty(a_z) && return (b, b_z)
    isempty(b_z) && return (a, a_z)

    la = length(a)
    lb = length(b)
    if a[end] == b[begin]
        total_len = la + lb - 1
        w = Vector{Variable}(undef, total_len)
        z = Vector{Int}(undef, total_len)
        @inbounds for k in 1:(la - 1)
            w[k] = a[k]
            z[k] = a_z[k]
        end

        w[la] = a[la]
        z[la] = a_z[la] + b_z[begin]

        @inbounds for k in 2:lb
            w[la + k - 1] = b[k]
            z[la + k - 1] = b_z[k]
        end
    else
        total_len = la + lb
        w = Vector{Variable}(undef, total_len)
        z = Vector{Int}(undef, total_len)
        @inbounds for k in 1:la
            w[k] = a[k]
            z[k] = a_z[k]
        end
        @inbounds for k in 1:lb
            w[la + k] = b[k]
            z[la + k] = b_z[k]
        end
    end
    return (w, z)
end

function Base.:(*)(x::Monomial, y::Monomial)
    isempty(x.z) && return y
    isempty(y.z) && return x
    return Monomial(_concat_var_expos(x.vars, x.z, y.vars, y.z)...)
end

"""
    neat_dot(x::Monomial, y::Monomial)

Computes the "neat dot" product of two monomials as star(x) * y.

# Arguments
- `x::Monomial`: First monomial
- `y::Monomial`: Second monomial

# Returns
- `Monomial`: Product of star(x) and y
"""
neat_dot(x::Monomial, y::Monomial) =
    Monomial(_concat_var_expos(reverse(x.vars), reverse(x.z), y.vars, y.z)...)

function _concat_var_expos3(
    a::Vector{Variable},
    a_z::Vector{Int},
    b::Vector{Variable},
    b_z::Vector{Int},
    c::Vector{Variable},
    c_z::Vector{Int},
)
    isempty(a_z) && return _concat_var_expos(b, b_z, c, c_z)
    isempty(b_z) && return _concat_var_expos(a, a_z, c, c_z)
    isempty(c_z) && return _concat_var_expos(a, a_z, b, b_z)

    la, lb, lc = length(a), length(b), length(c)

    if a[end] == b[begin] && b[end] == c[begin]
        total_len = la + lb + lc - 2
        w = Vector{Variable}(undef, total_len)
        z = Vector{Int}(undef, total_len)
        @inbounds for k in 1:(la - 1)
            w[k] = a[k]
            z[k] = a_z[k]
        end

        w[la] = a[la]
        z[la] = a_z[la] + b_z[begin]

        @inbounds for k in 2:(lb - 1)
            w[la + k - 1] = b[k]
            z[la + k - 1] = b_z[k]
        end

        if lb > 1
            w[la + lb - 1] = b[lb]
            z[la + lb - 1] = b_z[lb] + c_z[begin]
            @inbounds for k in 2:lc
                w[la + lb + k - 2] = c[k]
                z[la + lb + k - 2] = c_z[k]
            end
        else
            z[la] += c_z[begin]
            @inbounds for k in 2:lc
                w[la + k - 1] = c[k]
                z[la + k - 1] = c_z[k]
            end
        end

    elseif a[end] != b[begin] && b[end] == c[begin]
        total_len = la + lb + lc - 1
        w = Vector{Variable}(undef, total_len)
        z = Vector{Int}(undef, total_len)

        @inbounds for k in 1:la
            w[k] = a[k]
            z[k] = a_z[k]
        end

        @inbounds for k in 1:(lb - 1)
            w[la + k] = b[k]
            z[la + k] = b_z[k]
        end

        w[la + lb] = b[lb]
        z[la + lb] = b_z[lb] + c_z[begin]

        @inbounds for k in 2:lc
            w[la + lb + k - 1] = c[k]
            z[la + lb + k - 1] = c_z[k]
        end
    elseif a[end] == b[begin] && b[end] != c[begin]
        total_len = la + lb + lc - 1
        w = Vector{Variable}(undef, total_len)
        z = Vector{Int}(undef, total_len)

        @inbounds for k in 1:(la - 1)
            w[k] = a[k]
            z[k] = a_z[k]
        end

        w[la] = a[la]
        z[la] = a_z[la] + b_z[begin]

        @inbounds for k in 2:lb
            w[la + k - 1] = b[k]
            z[la + k - 1] = b_z[k]
        end

        @inbounds for k in 1:lc
            w[la + lb + k - 1] = c[k]
            z[la + lb + k - 1] = c_z[k]
        end
    else
        total_len = la + lb + lc
        w = Vector{Variable}(undef, total_len)
        z = Vector{Int}(undef, total_len)
        @inbounds for k in 1:la
            w[k] = a[k]
            z[k] = a_z[k]
        end
        @inbounds for k in 1:lb
            w[la + k] = b[k]
            z[la + k] = b_z[k]
        end
        @inbounds for k in 1:lc
            w[la + lb + k] = c[k]
            z[la + lb + k] = c_z[k]
        end
    end
    return (w, z)
end

function _neat_dot3(x::Monomial, y::Monomial, z::Monomial)
    return Monomial(
        _concat_var_expos3(reverse(x.vars), reverse(x.z), y.vars, y.z, z.vars, z.z)...
    )
end

Base.:(*)(a::Number, b::Polynomial) = Polynomial(a .* b.coeffs, b.monos)
Base.:(*)(a::Polynomial, b::Number) = Polynomial(b .* a.coeffs, a.monos)
Base.:(*)(a::Number, b::PolynomialLike) = a * Polynomial(b)
Base.:(*)(a::PolynomialLike, b::Number) = Polynomial(a) * b

Base.:(/)(a::PolynomialLike, b::Number) = Polynomial(a) * inv(b)
