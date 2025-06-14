
"""
    symmetric_canonicalize(mono::Monomial)

Canonicalizes a mono by taking the minimum between itself and its adjoint.

# Arguments
- `mono::Monomial`: The mono to canonicalize

# Returns
- `Monomial`: Canonicalized monomial (minimum of original and its star)
"""
function symmetric_canonicalize(mono::Monomial)
    isempty(mono.vars) && return mono
    return min(mono, star(mono))
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
    cyclic_canonicalize(mono::Monomial)

Canonicalizes a monomial using both cyclic and symmetric operations.
Finds the minimum among all cyclic shifts and their adjoints.
Chclic canonical is both cyclic and symmetric

# Arguments
- `mono::Monomial`: The mono to canonicalize

# Returns
- `Monomial`: Canonicalized monomial (minimum across all cyclic shifts and their stars)
"""
function cyclic_canonicalize(mono::Monomial)
    isempty(mono.vars) && return mono
    flatten_vars = mapreduce(
        idx -> fill(mono.vars[idx], mono.z[idx]), vcat, eachindex(mono.z)
    )
    flatten_z = ones(Int, sum(mono.z))
    return minimum(
        mapreduce(vcat, 1:sum(mono.z)) do shift
            shifted_mono = monomial(circshift!(flatten_vars, 1), circshift!(flatten_z, 1))
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
    sorted_unique(xs)

Returns sorted unique elements from a collection.

# Arguments
- `xs`: Collection to process

# Returns
- Sorted vector of unique elements
"""
sorted_unique(xs) = sort!(unique(xs))
sorted_unique!(xs) = sort!(unique!(xs))

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
    _comm(mono::Monomial, comm_gps::Vector{Vector{Variable}})

Projects a monomial onto commutative groups of variables while maintaining the
order of variables within each group.

We kept it as separate groups because symmetric canonicalize over the product of
groups vs product of symmetric canonicalize of each group is different


# Arguments
- `mono::Monomial`: The monomial to project
- `comm_gps::Vector{Vector{Variable}}`: Vector of sets defining commutative variable groups

# Returns
- `Vector{Monomial}`: Projections of the monomial onto each commutative group

# Example
```jldoctest; setup=:(using NCTSSoS.FastPolynomials; using NCTSSoS.FastPolynomials: _comm)
julia> @ncpolyvar x y; comm_gps = [[x], [y]]
2-element Vector{Vector{Variable}}:
 [x]
 [y]

julia> mono1 = x*y*x*y
x¹y¹x¹y¹

julia> _comm(mono1, comm_gps)
2-element Vector{Monomial}:
 x²
 y²
```
"""
function _comm(mono::Monomial, comm_gps::Vector{Vector{Variable}})
    map(comm_gps) do vars
        result = one(Monomial)
        for (var, expon) in zip(mono.vars, mono.z)
            if var in vars
                _mul_var!(result, var, expon, false)
            end
        end
        result
    end
end

# multiply a variable to a monomial
@inline function _mul_var!(result::Monomial, var::Variable, expo::Int, is_unipotent::Bool)
    # Q: do we need to consider commutative case?
    if is_unipotent
        iseven(expo) && return result
        if length(result.vars) == 0 || var != result.vars[end]   # new variable
            push!(result.vars, var)
            push!(result.z, expo)
        else
            pop!(result.vars)
            pop!(result.z)
        end
    else
        iszero(expo) && return result
        if length(result.vars) == 0 || var != result.vars[end]   # new variable
            push!(result.vars, var)
            push!(result.z, expo)
        else
            result.z[end] += expo
        end
    end
    return result
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
    result = one(Monomial)
    for (var, expo) in zip(mono.vars, mono.z)
        _mul_var!(result, var, expo, true)
    end
    return result
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
