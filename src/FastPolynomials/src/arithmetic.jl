const VarOrMono = Union{Variable,Monomial}

function Base.:(+)(a::VarOrMono, b::VarOrMono)
    return Polynomial(a) + Polynomial(b)
end

function Base.:(+)(a::VarOrMono, b::Polynomial{T}) where {T}
    return Polynomial(T, a) + b
end

function Base.:(+)(a::Polynomial{T}, b::VarOrMono) where {T}
    return a + Polynomial(T, b)
end

function Base.:(+)(a::Polynomial{T1}, b::Polynomial{T2}) where {T1<:Number,T2<:Number}
    T = promote_type(T1, T2)
    return Polynomial(T[a.coeffs; b.coeffs], [a.monos; b.monos])
end

function Base.:(-)(a::VarOrMono, b::VarOrMono)
    return Polynomial(a) - Polynomial(b)
end

function Base.:(-)(a::VarOrMono, b::Polynomial{T}) where {T<:Number}
    return Polynomial(T, a) - b
end

function Base.:(-)(a::Polynomial{T}, b::VarOrMono) where {T<:Number}
    return a - Polynomial(T, b)
end

function Base.:(-)(a::Polynomial{T1}, b::Polynomial{T2}) where {T1<:Number,T2<:Number}
    T = promote_type(T1, T2)
    return Polynomial(T[a.coeffs; -b.coeffs], [a.monos; b.monos])
end

function Base.:(-)(a::Polynomial{T1}, b::T2) where {T1,T2}
    return a - Polynomial(T1, b)
end

function Base.:(-)(a::Polynomial)
    return Polynomial(-a.coeffs, a.monos)
end

function Base.:(*)(a::Polynomial, b::Polynomial)
    return Polynomial(
        vec([ca * cb for (ca, cb) in Iterators.product(a.coeffs, b.coeffs)]),
        vec([ma * mb for (ma, mb) in Iterators.product(a.monos, b.monos)]),
    )
end

Base.:(*)(a::Polynomial{T}, b::VarOrMono) where {T} = a * Polynomial(T, b)
Base.:(*)(a::VarOrMono, b::Polynomial{T}) where {T} = Polynomial(T, a) * b

function Base.:(*)(a::T1, b::Polynomial) where {T1<:Number}
    return Polynomial(a .* b.coeffs, b.monos)
end

function Base.:(*)(a::T, b::VarOrMono) where {T<:Number}
    return a * Polynomial(T, b)
end

function ncmul(x::Monomial, y::Monomial)
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

function Base.:(*)(a::VarOrMono, b::VarOrMono)
    return ncmul(
        a isa Variable ? Monomial([a], [1]) : a, b isa Variable ? Monomial([b], [1]) : b
    )
end
