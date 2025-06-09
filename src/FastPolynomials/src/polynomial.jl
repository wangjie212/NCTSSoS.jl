"""
    Polynomial{T}

A polynomial represented as a sum of monomials with coefficients.
Maintains sorted, unique monomials with non-zero coefficients.

# Fields
- `coeffs::Vector{T}`: Non-zero coefficients
- `monos::Vector{Monomial}`: Sorted, unique monomials

# Invariants
- `monos` are always sorted and non-repeating
- `coeffs` are always non-zero
- Length of `coeffs` equals length of `monos`
"""
struct Polynomial{T}
    coeffs::Vector{T}
    # perhaps make it into an ordered set?
    monos::Vector{Monomial}

    function Polynomial(a::Vector{T}, x::Vector{Monomial}) where {T}
        length(a) == length(x) ||
            throw(ArgumentError("There should be as many coefficient than monomials"))
        sorted_x_idx = sortperm(x)
        sorted_x = x[sorted_x_idx]
        unique_monos = sort!(unique(x))
        unique_coeffs = map(unique_monos) do mono
            mapreduce(
                idx -> a[sorted_x_idx[idx]],
                +,
                searchsortedfirst(sorted_x, mono):searchsortedlast(sorted_x, mono);
                init=zero(T),
            )
        end
        nz_idx = findall(!iszero, unique_coeffs)
        return new{T}(unique_coeffs[nz_idx], unique_monos[nz_idx])
    end
end

"""
    Polynomial(a)

Creates a polynomial with Float64 coefficients.

# Arguments
- `a`: Variable, monomial, or scalar to convert to polynomial

# Returns
- `Polynomial{Float64}`: Polynomial with Float64 coefficients
"""
Polynomial(a) = Polynomial(Float64, a)

"""
    Polynomial(::Type{T}, a::Variable) where {T}

Creates a polynomial from a single variable with coefficient type T.

# Arguments
- `::Type{T}`: Coefficient type
- `a::Variable`: Variable to convert to polynomial

# Returns
- `Polynomial{T}`: Polynomial representing the variable with coefficient 1
"""
Polynomial(::Type{T}, a::Variable) where {T} = Polynomial([one(T)], [Monomial([a], [1])])

"""
    Polynomial(::Type{T}, a::Monomial) where {T}

Creates a polynomial from a single monomial with coefficient type T.

# Arguments
- `::Type{T}`: Coefficient type
- `a::Monomial`: Monomial to convert to polynomial

# Returns
- `Polynomial{T}`: Polynomial representing the monomial with coefficient 1
"""
Polynomial(::Type{T}, a::Monomial) where {T} = Polynomial([one(T)], [a])
"""
    Polynomial(::Type{T1}, a::T2) where {T1,T2}

Creates a constant polynomial from a scalar value.

# Arguments
- `::Type{T1}`: Target coefficient type
- `a::T2`: Scalar value

# Returns
- `Polynomial`: Constant polynomial with promoted coefficient type
"""
function Polynomial(::Type{T1}, a::T2) where {T1<:Number,T2<:Number}
    return Polynomial(promote_type(T1, T2), Monomial([], []))
end

Polynomial(a::Polynomial) = a

function Base.convert(::Type{Polynomial{T}}, a::Variable) where {T}
    return Polynomial([one(T)], [Monomial([a], [one(T)])])
end
Base.convert(::Type{Polynomial{T}}, a::Monomial) where {T} = Polynomial([one(T)], [a])

"""
    variables(p::Polynomial)

Extracts all variables appearing in a polynomial.

# Arguments
- `p::Polynomial`: The polynomial

# Returns
- Set of all variables in the polynomial
"""
variables(p::Polynomial) = union(variables.(p.monos)...)

function Base.show(io::IO, obj::Polynomial)
    return print_object(io, obj; multiline=true)
end

function Base.show(io::IO, mime::MIME"text/plain", obj::Polynomial)
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

"""
    print_object(io, obj; multiline)

Prints a Polynomial object to an IO stream.

# Arguments
- `io::IO`: The output stream
- `obj::Polynomial`: The polynomial to print
- `multiline::Bool`: Whether to use expanded format

# Returns
- Nothing (prints to IO stream)
"""
function print_object(io::IO, obj::Polynomial; multiline::Bool)
    if multiline
        return print(io, join(
            map(zip(obj.coeffs, obj.monos)) do (coef, mono)
                "$(coef) * $(string(mono))"
            end,
            " + ",
        ))
    else
        return Base.show_default(io, obj)
    end
end

"""
    Base.hash(p::Polynomial, u::UInt)

Computes hash value for a polynomial based on its coefficients and monomials.

# Arguments
- `p::Polynomial`: The polynomial to hash
- `u::UInt`: Hash seed value

# Returns
- `UInt`: Hash value combining coefficients and monomials
"""
function Base.hash(p::Polynomial, u::UInt)
    return hash(p.coeffs, hash(p.monos, u))
end

maxdegree(p::Polynomial) = maximum(degree.(p.monos))

function Base.:(^)(p::Polynomial, n::Int)
    n == 0 && return one(typeof(p))
    n == 1 && return p
    return p * (p^(n - 1))
end

Base.one(::Polynomial{T}) where {T} = Polynomial([one(T)], [Monomial([], [])])
Base.zero(::Polynomial{T}) where {T} = Polynomial([zero(T)], [Monomial([], [])])

coefficients(p::Polynomial) = p.coeffs
monomials(p::Polynomial) = p.monos

terms(p::Polynomial) = zip(p.coeffs, p.monos)
