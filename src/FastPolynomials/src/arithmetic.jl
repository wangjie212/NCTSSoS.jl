"""
    VarOrMono

Type alias for objects that can be either a Variable or a Monomial.
Used in arithmetic operations to allow flexible input types.
"""
const VarOrMono = Union{Variable,Monomial}

"""
    Base.:(+)(a::VarOrMono, b::VarOrMono)

Adds two variables/monomials by converting them to polynomials.

# Arguments
- `a::VarOrMono`: First variable or monomial
- `b::VarOrMono`: Second variable or monomial

# Returns
- `Polynomial`: Sum of the two inputs as a polynomial
"""
function Base.:(+)(a::VarOrMono, b::VarOrMono)
    return Polynomial(a) + Polynomial(b)
end

"""
    Base.:(+)(a::VarOrMono, b::Polynomial{T}) where {T}

Adds a variable/monomial to a polynomial.

# Arguments
- `a::VarOrMono`: Variable or monomial to add
- `b::Polynomial{T}`: Polynomial to add to

# Returns
- `Polynomial{T}`: Sum with promoted coefficient type
"""
function Base.:(+)(a::VarOrMono, b::Polynomial{T}) where {T}
    return Polynomial(T, a) + b
end

"""
    Base.:(+)(a::Polynomial{T}, b::VarOrMono) where {T}

Adds a polynomial to a variable/monomial.

# Arguments
- `a::Polynomial{T}`: Polynomial to add
- `b::VarOrMono`: Variable or monomial to add

# Returns
- `Polynomial{T}`: Sum with coefficient type T
"""
function Base.:(+)(a::Polynomial{T}, b::VarOrMono) where {T}
    return a + Polynomial(T, b)
end

"""
    Base.:(+)(a::Polynomial{T1}, b::Polynomial{T2}) where {T1<:Number,T2<:Number}

Adds two polynomials with potentially different coefficient types.

# Arguments
- `a::Polynomial{T1}`: First polynomial
- `b::Polynomial{T2}`: Second polynomial

# Returns
- `Polynomial`: Sum with promoted coefficient type combining all terms
"""
function Base.:(+)(a::Polynomial{T1}, b::Polynomial{T2}) where {T1<:Number,T2<:Number}
    T = promote_type(T1, T2)
    return Polynomial(T[a.coeffs; b.coeffs], [a.monos; b.monos])
end

"""
    Base.:(-)(a::VarOrMono, b::VarOrMono)

Subtracts two variables/monomials by converting them to polynomials.

# Arguments
- `a::VarOrMono`: Variable or monomial (minuend)
- `b::VarOrMono`: Variable or monomial (subtrahend)

# Returns
- `Polynomial`: Difference of the two inputs as a polynomial
"""
function Base.:(-)(a::VarOrMono, b::VarOrMono)
    return Polynomial(a) - Polynomial(b)
end

"""
    Base.:(-)(a::VarOrMono, b::Polynomial{T}) where {T<:Number}

Subtracts a polynomial from a variable/monomial.

# Arguments
- `a::VarOrMono`: Variable or monomial (minuend)
- `b::Polynomial{T}`: Polynomial (subtrahend)

# Returns
- `Polynomial{T}`: Difference with coefficient type T
"""
function Base.:(-)(a::VarOrMono, b::Polynomial{T}) where {T<:Number}
    return Polynomial(T, a) - b
end

"""
    Base.:(-)(a::Polynomial{T}, b::VarOrMono) where {T<:Number}

Subtracts a variable/monomial from a polynomial.

# Arguments
- `a::Polynomial{T}`: Polynomial (minuend)
- `b::VarOrMono`: Variable or monomial (subtrahend)

# Returns
- `Polynomial{T}`: Difference with coefficient type T
"""
function Base.:(-)(a::Polynomial{T}, b::VarOrMono) where {T<:Number}
    return a - Polynomial(T, b)
end

"""
    Base.:(-)(a::Polynomial{T1}, b::Polynomial{T2}) where {T1<:Number,T2<:Number}

Subtracts two polynomials with potentially different coefficient types.

# Arguments
- `a::Polynomial{T1}`: First polynomial (minuend)
- `b::Polynomial{T2}`: Second polynomial (subtrahend)

# Returns
- `Polynomial`: Difference with promoted coefficient type, negating coefficients of b
"""
function Base.:(-)(a::Polynomial{T1}, b::Polynomial{T2}) where {T1<:Number,T2<:Number}
    T = promote_type(T1, T2)
    return Polynomial(T[a.coeffs; -b.coeffs], [a.monos; b.monos])
end

"""
    Base.:(-)(a::Polynomial{T1}, b::T2) where {T1,T2}

Subtracts a scalar from a polynomial.

# Arguments
- `a::Polynomial{T1}`: Polynomial (minuend)
- `b::T2`: Scalar value (subtrahend)

# Returns
- `Polynomial`: Difference treating scalar as constant polynomial
"""
function Base.:(-)(a::Polynomial{T1}, b::T2) where {T1,T2<:Number}
    return a - Polynomial(T1, b)
end

"""
    Base.:(-)(a::Polynomial)

Negates a polynomial by negating all coefficients.

# Arguments
- `a::Polynomial`: Polynomial to negate

# Returns
- `Polynomial`: Negated polynomial with same monomials but negated coefficients
"""
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

"""
    Base.:(*)(a::Polynomial{T}, b::VarOrMono) where {T}

Multiplies a polynomial by a variable/monomial.

# Arguments
- `a::Polynomial{T}`: Polynomial
- `b::VarOrMono`: Variable or monomial

# Returns
- `Polynomial{T}`: Product with coefficient type T
"""
Base.:(*)(a::Polynomial{T}, b::VarOrMono) where {T} = a * Polynomial(T, b)

"""
    Base.:(*)(a::VarOrMono, b::Polynomial{T}) where {T}

Multiplies a variable/monomial by a polynomial.

# Arguments
- `a::VarOrMono`: Variable or monomial
- `b::Polynomial{T}`: Polynomial

# Returns
- `Polynomial{T}`: Product with coefficient type T
"""
Base.:(*)(a::VarOrMono, b::Polynomial{T}) where {T} = Polynomial(T, a) * b

"""
    Base.:(*)(a::T1, b::Polynomial) where {T1<:Number}

Multiplies a polynomial by a scalar coefficient.

# Arguments
- `a::T1`: Scalar multiplier
- `b::Polynomial`: Polynomial to multiply

# Returns
- `Polynomial`: Polynomial with all coefficients multiplied by the scalar
"""
function Base.:(*)(a::T1, b::Polynomial) where {T1<:Number}
    return Polynomial(a .* b.coeffs, b.monos)
end

"""
    Base.:(*)(a::T, b::VarOrMono) where {T<:Number}

Multiplies a variable/monomial by a scalar coefficient.

# Arguments
- `a::T`: Scalar multiplier
- `b::VarOrMono`: Variable or monomial

# Returns
- `Polynomial{T}`: Polynomial with the variable/monomial scaled by the coefficient
"""
function Base.:(*)(a::T, b::VarOrMono) where {T<:Number}
    return a * Polynomial(T, b)
end

"""
    ncmul(x::Monomial, y::Monomial)

Non-commutative multiplication of two monomials.
Concatenates monomials, combining exponents when adjacent variables match.

# Arguments
- `x::Monomial`: First monomial (left factor)
- `y::Monomial`: Second monomial (right factor)

# Returns
- `Monomial`: Product monomial respecting non-commutative variable order
"""
function ncmul(x::Monomial, y::Monomial)
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

"""
    Base.:(*)(a::VarOrMono, b::VarOrMono)

Non-commutative multiplication of variables/monomials.

# Arguments
- `a::VarOrMono`: First variable or monomial
- `b::VarOrMono`: Second variable or monomial

# Returns
- `Monomial`: Product using non-commutative multiplication rules
"""
function Base.:(*)(a::VarOrMono, b::VarOrMono)
    return ncmul(
        a isa Variable ? Monomial([a], [1]) : a, b isa Variable ? Monomial([b], [1]) : b
    )
end
