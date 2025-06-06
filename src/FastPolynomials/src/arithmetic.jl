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

function Base.:(+)(a::Polynomial{T1}, b::Polynomial{T2}) where {T1<:Number,T2<:Number}
    T = promote_type(T1, T2)
    return Polynomial(T[a.coeffs; b.coeffs], [a.monos; b.monos])
end

function Base.:(+)(a::Polynomial, b::Monomial)
    return Polynomial([a.coeffs; one(eltype(a.coeffs))], [a.monos; b])
end

function Base.:(*)(a::Number, b::Polynomial{T}) where {T}
    return Polynomial(a .* b.coeffs, b.monos)
end

function Base.:(*)(a::T, b::Variable) where {T<:Number}
    return Polynomial([a], [Monomial([b], [1])])
end
