struct Polynomial{T}
	coeffs::Vector{T}
	monos::Vector{Monomial}

    function Polynomial(
        a::Vector{T},
        x::Vector{Monomial},
    ) where T 
        length(a) == length(x) || throw(
            ArgumentError("There should be as many coefficient than monomials"),
        )
		nz_idx = findall(!iszero, a)
		return new{T}(a[nz_idx], x[nz_idx])
    end
end

