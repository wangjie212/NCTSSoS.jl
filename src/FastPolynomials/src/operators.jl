"""
    Base.:(^)(a::Variable, expo::Int)

Raises a variable to a non-negative integer power, creating a monomial.

# Arguments
- `a::Variable`: The variable to raise to a power
- `expo::Int`: Non-negative integer exponent

# Returns
- `Monomial`: Empty monomial if expo is 0, otherwise monomial with variable raised to expo

# Throws
- `AssertionError`: If expo is negative
"""
function Base.:(^)(a::Variable, expo::Int)
    @assert expo >= 0 "Exponent must be non-negative."
    return iszero(expo) ? Monomial([], []) : Monomial([a], [expo])
end
