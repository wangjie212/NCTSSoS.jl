function Base.:(^)(a::Variable, expo::Int)
    @assert expo >= 0 "Exponent must be non-negative."
    return iszero(expo) ? Monomial([], []) : Monomial([a], [expo])
end

function Base.:(*)(coef::T, a::Monomial) where {T<:Number}
    return Polynomial([coef], [a])
end

function Base.:(+)(a::Monomial, b::Monomial)
    return Polynomial([1.0, 1.0], [a, b])
end

function Base.:(*)(a::Variable, b::Variable)
    return Monomial([a, b], [1, 1])
end

function Base.:(*)(a::Monomial, b::Variable)
    return Monomial([a.vars; b], [a.z; 1])
end

function Base.:(*)(a::Variable, b::Monomial)
    return Monomial([a; b.vars], [1; b.z])
end
