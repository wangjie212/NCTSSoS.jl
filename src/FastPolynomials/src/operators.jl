function Base.:(^)(a::Variable, expo::Int)
    @assert expo >= 0 "Exponent must be non-negative."
    return iszero(expo) ? Monomial([], []) : Monomial([a], [expo])
end
