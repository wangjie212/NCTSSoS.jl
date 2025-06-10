"""
    star(m::Monomial)

Computes the adjoint (star) of a monomial by reversing variable order and exponents.

# Arguments
- `m::Monomial`: The monomial to compute the adjoint of

# Returns
- `Monomial`: Adjoint monomial with reversed variables and exponents
"""
function star(m::Monomial)
    return Monomial(reverse(m.vars), reverse(m.z))
end

"""
    symmetric_canonicalize(monomial::Monomial)

Canonicalizes a monomial by taking the minimum between itself and its adjoint.

# Arguments
- `monomial::Monomial`: The monomial to canonicalize

# Returns
- `Monomial`: Canonicalized monomial (minimum of original and its star)
"""
function symmetric_canonicalize(monomial::Monomial)
    isempty(monomial.vars) && return monomial
    return min(monomial, star(monomial))
end

"""
    symmetric_canonicalize(poly::Polynomial)

Canonicalizes a polynomial by applying symmetric canonicalization to all monomials.

# Arguments
- `poly::Polynomial`: The polynomial to canonicalize

# Returns
- `Polynomial`: Canonicalized polynomial with conjugated coefficients and canonicalized monomials
"""
function symmetric_canonicalize(poly::Polynomial)
    return Polynomial(conj.(poly.coeffs), symmetric_canonicalize.(poly.monos))
end

"""
    cyclic_canonicalize(monomial::Monomial)

Canonicalizes a monomial using both cyclic and symmetric operations.
Finds the minimum among all cyclic shifts and their adjoints.
Chclic canonical is both cyclic and symmetric

# Arguments
- `monomial::Monomial`: The monomial to canonicalize

# Returns
- `Monomial`: Canonicalized monomial (minimum across all cyclic shifts and their stars)
"""
function cyclic_canonicalize(monomial::Monomial)
    isempty(monomial.vars) && return monomial
    flatten_vars = mapreduce(
        idx -> fill(monomial.vars[idx], monomial.z[idx]), vcat, eachindex(monomial.z)
    )
    flatten_z = ones(Int, sum(monomial.z))
    return minimum(
        mapreduce(vcat, 1:sum(monomial.z)) do shift
            shifted_mono = Monomial(circshift!(flatten_vars, 1), circshift!(flatten_z, 1))
            [shifted_mono, star(shifted_mono)]
        end,
    )
end

"""
    cyclic_canonicalize(poly::Polynomial)

Canonicalizes a polynomial by applying cyclic canonicalization to all monomials.

# Arguments
- `poly::Polynomial`: The polynomial to canonicalize

# Returns
- `Polynomial`: Canonicalized polynomial with same coefficients and canonicalized monomials
"""
function cyclic_canonicalize(poly::Polynomial)
    return Polynomial(poly.coeffs, cyclic_canonicalize.(poly.monos))
end



"""
    support(poly::Polynomial{T}, canonicalize::Function) where {T}

Computes the support of a polynomial after canonicalization.

# Arguments
- `poly::Polynomial{T}`: The polynomial
- `canonicalize::Function`: Function to canonicalize support

# Returns
- `Vector{Monomial}`: Unique canonicalized monomials from the polynomial
"""
function support(poly::Polynomial{T}, canonicalize::Function) where {T}
    return unique!(canonicalize.(poly.monos))
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
function neat_dot(x::Monomial, y::Monomial)
    return star(x) * y
end


"""
    sorted_unique(xs)

Returns sorted unique elements from a collection.

# Arguments
- `xs`: Collection to process

# Returns
- Sorted vector of unique elements
"""
sorted_unique(xs) = sort(unique(xs))

"""
    sorted_union(xs...)

Returns sorted union of multiple collections.

# Arguments
- `xs...`: Variable number of collections

# Returns
- Sorted vector containing union of all input collections
"""
sorted_union(xs...) = sort(union(xs...))

"""
    _comm(mono::Monomial, comm_gps::Vector{Set{Variable}})

Projects a monomial onto commutative groups of variables while maintaining the
order of variables within each group.

# Arguments
- `mono::Monomial`: The monomial to project
- `comm_gps::Vector{Set{Variable}}`: Vector of sets defining commutative variable groups

# Returns
- `Vector{Monomial}`: Projections of the monomial onto each commutative group

# Example
```jldoctest; setup=:(using NCTSSoS.FastPolynomials; using NCTSSoS.FastPolynomials: _comm)
julia> @ncpolyvar x y; comm_gps = [Set([x]),Set([y])]
2-element Vector{Set{Variable}}:
 Set([x])
 Set([y])

julia> mono1 = x*y*x*y
x¹y¹x¹y¹

julia> _comm(mono1, comm_gps)
2-element Vector{Monomial}:
 x²
 y²
```
"""
function _comm(mono::Monomial, comm_gps::Vector{Set{Variable}})
    map(comm_gps) do vars
        prod(zip(mono.vars, mono.z); init=Monomial([], [])) do (var, expon)
            var in vars ? var^expon : var^(zero(expon))
        end
    end
end

"""
    _unipotent(mono::Monomial)

Applies unipotent transformation to a monomial by reducing exponents modulo 2 iteratively.

# Arguments
- `mono::Monomial`: The monomial to transform

# Returns
- `Monomial`: Unipotent form of the monomial with all exponents reduced to 0 or 1
"""
function _unipotent(mono::Monomial)
    isempty(mono.vars) && return mono
    prev_mono = mono
    local cur_mono
    while true
        cur_mono = prod(zip(prev_mono.vars, prev_mono.z); init=one(Monomial)) do (var, expo)
            var^(expo % 2)
        end
        cur_mono == prev_mono && break
        prev_mono = cur_mono
    end
    return cur_mono
end


"""
    _projective(mono::Monomial)

Applies projective transformation to a monomial by setting all non-zero exponents to 1.

# Arguments
- `mono::Monomial`: The monomial to transform

# Returns
- `Monomial`: Projective form with all non-zero exponents set to 1
"""
function _projective(mono::Monomial)
    prod(zip(mono.vars, mono.z); init=one(Monomial)) do (var, expo)
        var^(iszero(expo) ? expo : one(expo))
    end
end
