Base.cmp(a::Variable, b::Variable) = a.name == b.name ? 0 : (a.name < b.name ? -1 : 1)

Base.isless(a::Variable, b::Variable) = cmp(a, b) < 0

function Base.in(a::Variable, collection::Vector{Variable})
    return searchsortedfirst(collection, a; lt=cmp) != 0
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

Base.isless(a::Monomial, b::Monomial) = cmp(a, b) == -1

Base.:(==)(a::Monomial, b::Monomial) = iszero(cmp(a, b))

# Comparison of Term
function (==)(p::Polynomial{T}, q::Polynomial{T}) where {T}
    length(p.monos) != length(q.monos) && return false
    map(zip(p.monos, q.monos)) do (mono1, mono2)
        mono1 != mono2 && return false
    end
    map(zip(p.coeffs, q.coeffs)) do (mono1, mono2)
        mono1 != mono2 && return false
    end
    return true
end

# FIXME
function Base.isapprox(
    p::Polynomial{S},
    q::Polynomial{T};
    rtol::Real=Base.rtoldefault(S, T, 0),
    atol::Real=0.0,
    ztol::Real=iszero(atol) ? Base.rtoldefault(S, T, 0) : atol,
) where {S,T}
    i = j = 1
    while i <= length(p.monos) || j <= length(q.monos)
        if i > length(p.x) || (j <= length(q.x) && q.x[j] < p.x[i])
            if !isapproxzero(q.a[j]; ztol=ztol)
                return false
            end
            j += 1
        elseif j > length(q.x) || p.x[i] < q.x[j]
            if !isapproxzero(p.a[i]; ztol=ztol)
                return false
            end
            i += 1
        else
            if !isapprox(p.a[i], q.a[j]; rtol=rtol, atol=atol)
                return false
            end
            i += 1
            j += 1
        end
    end
    return true
end
