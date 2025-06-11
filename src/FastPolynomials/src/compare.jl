Base.cmp(a::Variable, b::Variable) = (a.name == b.name) ? 0 : (a.name < b.name ? -1 : 1)

Base.isless(a::Variable, b::Variable) = cmp(a, b) < 0

function Base.in(a::Variable, collection::Vector{Variable})
    if !issorted(collection)
        sort!(collection)
    end
    return !isempty(searchsorted(collection, a))
end

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

function binary_search(collection::Vector{Monomial}, target::Monomial)
    left, right = 1, length(collection)
    while left < right
        mid = (left + right) ÷ 2
        cmp_result = cmp(target, collection[mid])
        if cmp_result < 0
            right = mid
        elseif cmp_result > 0
            left = mid + 1
        else
            return mid
        end
    end
    return 0
end

function Base.in(a::Monomial, collection::Vector{Monomial})
    if !issorted(collection)
        sort!(collection)
    end
    return !iszero(binary_search(collection, a))
end

Base.isless(a::Monomial, b::Monomial) = cmp(a, b) == -1

Base.:(==)(a::Monomial, b::Monomial) = iszero(cmp(a, b))

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

function Base.isapprox(p::Polynomial{S}, q::Polynomial{T}; atol::Real=0.0) where {S,T}
    p_nz_idcs = findall(x -> !isapprox(0.0, p.coeffs[x]; atol=atol), 1:length(p.coeffs))
    q_nz_idcs = findall(x -> !isapprox(0.0, q.coeffs[x]; atol=atol), 1:length(q.coeffs))

    length(p_nz_idcs) != length(q_nz_idcs) && return false

    p.monos[p_nz_idcs] != q.monos[q_nz_idcs] && return false
    isapprox(p.coeffs[p_nz_idcs], q.coeffs[q_nz_idcs]; atol=atol) || return false
    return true
end
