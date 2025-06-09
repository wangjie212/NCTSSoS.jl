"""
    Base.cmp(a::Variable, b::Variable)

Compares two variables lexicographically by their names.

# Arguments
- `a::Variable`: First variable to compare
- `b::Variable`: Second variable to compare

# Returns
- `Int`: -1 if a < b, 0 if a == b, 1 if a > b
"""
Base.cmp(a::Variable, b::Variable) = (a.name == b.name) ? 0 : (a.name < b.name ? -1 : 1)

"""
    Base.isless(a::Variable, b::Variable)

Determines if one variable is lexicographically less than another.

# Arguments
- `a::Variable`: First variable
- `b::Variable`: Second variable

# Returns
- `Bool`: True if a's name is lexicographically less than b's name
"""
Base.isless(a::Variable, b::Variable) = cmp(a, b) < 0

"""
    Base.in(a::Variable, collection::Vector{Variable})

Checks if a variable is present in a sorted collection of variables.

# Arguments
- `a::Variable`: Variable to search for
- `collection::Vector{Variable}`: Sorted vector of variables

# Returns
- `Bool`: True if variable is found in the collection
"""
function Base.in(a::Variable, collection::Vector{Variable})
    @assert issorted(collection) "Collection is unsorted"
    return searchsortedfirst(collection, a) != 0
end

"""
    Base.cmp(a::Monomial, b::Monomial)

Compares two monomials using graded lexicographic ordering.
First compares by total degree, then lexicographically by variables and exponents.

# Arguments
- `a::Monomial`: First monomial to compare
- `b::Monomial`: Second monomial to compare

# Returns
- `Int`: -1 if a < b, 0 if a == b, 1 if a > b
"""
function Base.cmp(a::Monomial, b::Monomial)
    degree(a) != degree(b) && return degree(a) < degree(b) ? -1 : 1
    a_idx, b_idx = 1, 1
    while a_idx <= length(a.vars) && b_idx <= length(b.vars)
        iszero(a.z[a_idx]) && (a_idx += 1; continue)
        iszero(b.z[b_idx]) && (b_idx += 1; continue)
        a.vars[a_idx] != b.vars[b_idx] && return cmp(a.vars[a_idx], b.vars[b_idx])
        a.z[a_idx] != b.z[b_idx] && return cmp(a.z[a_idx], b.z[b_idx])
        a_idx += 1
        b_idx += 1
    end
    return 0
end

"""
    Base.isless(a::Monomial, b::Monomial)

Determines if one monomial is less than another using graded lexicographic ordering.

# Arguments
- `a::Monomial`: First monomial
- `b::Monomial`: Second monomial

# Returns
- `Bool`: True if a is less than b in graded lexicographic order
"""
Base.isless(a::Monomial, b::Monomial) = cmp(a, b) == -1

"""
    Base.:(==)(a::Monomial, b::Monomial)

Tests equality between two monomials.

# Arguments
- `a::Monomial`: First monomial
- `b::Monomial`: Second monomial

# Returns
- `Bool`: True if monomials are equal (same variables and exponents)
"""
Base.:(==)(a::Monomial, b::Monomial) = iszero(cmp(a, b))

"""
    Base.:(==)(p::Polynomial{T}, q::Polynomial{T}) where {T}

Tests exact equality between two polynomials of the same coefficient type.

# Arguments
- `p::Polynomial{T}`: First polynomial
- `q::Polynomial{T}`: Second polynomial

# Returns
- `Bool`: True if polynomials have identical monomials and coefficients
"""
function Base.:(==)(p::Polynomial{T}, q::Polynomial{T}) where {T}
    length(p.monos) != length(q.monos) && return false
    for (mono1, mono2) in zip(p.monos, q.monos)
        mono1 != mono2 && return false
    end
    for (coef1, coef2) in zip(p.coeffs, q.coeffs)
        coef1 != coef2 && return false
    end
    return true
end

"""
    Base.isapprox(p::Polynomial{S}, q::Polynomial{T}; atol::Real=0.0) where {S,T}

Tests approximate equality between two polynomials, ignoring near-zero coefficients.

# Arguments
- `p::Polynomial{S}`: First polynomial
- `q::Polynomial{T}`: Second polynomial
- `atol::Real`: Absolute tolerance for coefficient comparison (default: 0.0)

# Returns
- `Bool`: True if polynomials are approximately equal after removing near-zero terms
"""
function Base.isapprox(p::Polynomial{S}, q::Polynomial{T}; atol::Real=0.0) where {S,T}
    p_nz_idcs = findall(x -> !isapprox(0.0, p.coeffs[x]; atol=atol), 1:length(p.coeffs))
    q_nz_idcs = findall(x -> !isapprox(0.0, q.coeffs[x]; atol=atol), 1:length(q.coeffs))

    length(p_nz_idcs) != length(q_nz_idcs) && return false

    p.monos[p_nz_idcs] != q.monos[q_nz_idcs] && return false
    isapprox(p.coeffs[p_nz_idcs], q.coeffs[q_nz_idcs]; atol=atol) || return false
    return true
end
