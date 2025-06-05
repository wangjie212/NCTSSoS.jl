# WARNING:
# ALWAYS needs to guarantee `monos` are sorted
# ALWAYS needs to guarantee `coeffs` are non-zero
struct Polynomial{T}
    coeffs::Vector{T}
    # perhaps make it into an ordered set?
    monos::Vector{Monomial}

    function Polynomial(a::Vector{T}, x::Vector{Monomial}) where {T}
        length(a) == length(x) ||
            throw(ArgumentError("There should be as many coefficient than monomials"))
        nz_idx = findall(!iszero, a)
        sort!(nz_idx; by=idx -> x[idx])
        return new{T}(a[nz_idx], x[nz_idx])
    end
end
