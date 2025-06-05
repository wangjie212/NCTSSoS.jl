# I forbid the multiplication of two variables
# Base.:(*)(a::Variable, b::variable)

function Base.:(*)(x::Monomial, y::Monomial)
    i = findlast(z -> z > 0, x.z)
    isnothing(i) && return y
    j = findfirst(z -> z > 0, y.z)
    isnothing(j) && return x

    if x.vars[i] == y.vars[j]
        w = [x.vars[1:i]; y.vars[(j + 1):end]]
        z = [x.z[1:(i - 1)]; x.z[i] + y.z[j]; y.z[(j + 1):end]]
    else
        w = [x.vars[1:i]; y.vars[j:end]]
        z = [x.z[1:i]; y.z[j:end]]
    end
    return Monomial(w, z)
end
