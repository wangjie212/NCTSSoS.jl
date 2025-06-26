Base.:(*)(coef::Number, a::StateWord) = StatePolynomial([coef], [a])
Base.:(*)(coef::Number, a::NCStateWord) = NCStatePolynomial([coef], [a])

"""
    StatePolynomial{T}

A polynomial in state words with coefficients of type T.
Automatically combines like terms and maintains sorted unique state words.

# Fields
- `coeffs::Vector{T}`: Coefficients for each state word
- `state_words::Vector{StateWord}`: Sorted unique state words

# Constructor
Creates a StatePolynomial by combining coefficients for identical state words.
"""
struct StatePolynomial{T,ST} <: AbstractPolynomial{T}
    coeffs::Vector{T}
    state_words::Vector{StateWord{ST}}
    function StatePolynomial(coeffs::Vector{T}, state_words::Vector{StateWord{ST}}) where {T,ST}
        uniq_state_words = sorted_unique(state_words)
        uniq_coeffs = zeros(T, length(uniq_state_words))
        for (coef, sw) in zip(coeffs, state_words)
            idx = searchsortedfirst(uniq_state_words, sw)
            uniq_coeffs[idx] += coef
        end
        nz_idcs = filter(a -> !iszero(uniq_coeffs[a]), 1:length(uniq_coeffs))
        return new{T,ST}(uniq_coeffs[nz_idcs], uniq_state_words[nz_idcs])
    end
end

variables(sp::StatePolynomial) = sorted_union(variables.(sp.state_words)...)

degree(sp::StatePolynomial) = mapreduce(degree, max, sp.state_words)

monomials(sp::StatePolynomial) = sp.state_words

coefficients(sp::StatePolynomial) = sp.coeffs

function Base.show(io::IO, obj::StatePolynomial)
    return print_object(io, obj; multiline=false)
end

function Base.show(io::IO, mime::MIME"text/plain", obj::StatePolynomial)
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

function Base.string(obj::StatePolynomial)
    return join(
        map(zip(obj.coeffs, obj.state_words)) do (coef, sw)
            string(coef) * " * " * string(sw)
        end,
        " + ",
    )
end

function print_object(io::IO, obj::StatePolynomial; multiline::Bool)
    return multiline ? print(io, string(obj)) : Base.show_default(io, obj)
end

function Base.:(==)(a::StatePolynomial, b::StatePolynomial)
    a.state_words != b.state_words && return false
    a.coeffs != b.coeffs && return false
    return true
end

Base.hash(a::StatePolynomial, u::UInt) = hash(hash.(a.coeffs, u), hash.(a.state_words, u))

function Base.:(*)(a::StatePolynomial{T}, b::StatePolynomial{T}) where {T}
    return StatePolynomial(
        vec([ac * bc for (ac, bc) in Iterators.product(a.coeffs, b.coeffs)]),
        vec([asw * bsw for (asw, bsw) in Iterators.product(a.state_words, b.state_words)]),
    )
end

function Base.:(*)(a::StatePolynomial{T,ST}, b::Monomial) where {T,ST}
    return NCStatePolynomial(a.coeffs, [sw * b for sw in a.state_words])
end

function Base.:(*)(n::Number, a::StatePolynomial{T}) where {T}
    return StatePolynomial(T(n) .* a.coeffs, a.state_words)
end

function Base.:(*)(a::StatePolynomial, b::StateWord)
    return StatePolynomial(a.coeffs, [sw * b for sw in a.state_words])
end

function Base.:(+)(a::StatePolynomial{T1}, b::StatePolynomial{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    return StatePolynomial(T[a.coeffs; b.coeffs], [a.state_words; b.state_words])
end

function Base.:(+)(a::StatePolynomial{T}, b::StateWord) where {T}
    return StatePolynomial([a.coeffs; one(T)], [a.state_words; b])
end

function Base.:(+)(a::StateWord, b::StateWord)
    return StatePolynomial([one(Float64); one(Float64)], [a; b])
end

function Base.:(+)(a::Number, b::StateWord)
    return StatePolynomial([a; one(a)], [one(StateWord); b])
end

function Base.:(-)(a::StateWord, b::StateWord)
    return StatePolynomial([one(Float64); one(Float64)], [a; b])
end

function Base.:(-)(a::StatePolynomial{T}, b::StateWord) where {T}
    return StatePolynomial([a.coeffs; -one(T)], [a.state_words; b])
end

function Base.:(-)(a::StatePolynomial, b::StatePolynomial)
    return StatePolynomial([a.coeffs; -b.coeffs], [a.state_words; b.state_words])
end

Base.one(::StatePolynomial{T}) where {T} = StatePolynomial([one(T)], [one(StateWord)])
Base.one(::Type{StatePolynomial{T}}) where {T} = StatePolynomial([one(T)], [one(StateWord)])

function Base.zero(::StatePolynomial{T}) where {T}
    return StatePolynomial([zero(T)], [one(StateWord)])
end

terms(sp::StatePolynomial) = zip(sp.coeffs, sp.state_words)

"""
    NCStatePolynomial{T}

A polynomial in non-commutative state words with coefficients of type T.
Automatically combines like terms and maintains sorted unique NC state words.

# Fields
- `coeffs::Vector{T}`: Coefficients for each NC state word
- `nc_state_words::Vector{NCStateWord}`: Sorted unique non-commutative state words

# Constructor
Creates an NCStatePolynomial by combining coefficients for identical NC state words.
"""
# T: type of coefficient
struct NCStatePolynomial{T,ST} <: AbstractPolynomial{T}
    coeffs::Vector{T}
    nc_state_words::Vector{NCStateWord{ST}}
    function NCStatePolynomial(
        coeffs::Vector{T}, nc_state_words::Vector{NCStateWord{ST}}
    ) where {T,ST}
        @assert length(coeffs) == length(nc_state_words) "length of coeffs and nc_state_words must be the same, got $(length(coeffs)) and $(length(nc_state_words))"
        @assert issorted(nc_state_words) "nc_state_words must be sorted, got: $(nc_state_words)"
        @assert allunique(nc_state_words) "nc_state_words must be unique, got: $(nc_state_words)"
        return new{T,ST}(coeffs, nc_state_words)
    end
end
# the safer way to construct an NCStatePolynomial
function ncstatepoly(coeffs::Vector{T}, nc_state_words::Vector{NCStateWord{ST}}) where {T,ST}
    uniq_nc_state_words = sorted_unique(nc_state_words)
    uniq_coeffs = zeros(T, length(uniq_nc_state_words))
    for (coef, sw) in zip(coeffs, nc_state_words)
        idx = searchsortedfirst(uniq_nc_state_words, sw)
        uniq_coeffs[idx] += coef
    end
    nz_idcs = filter(a -> !iszero(uniq_coeffs[a]), 1:length(uniq_coeffs))
    return NCStatePolynomial(uniq_coeffs[nz_idcs], uniq_nc_state_words[nz_idcs])
end


function Base.show(io::IO, obj::NCStatePolynomial)
    return print_object(io, obj; multiline=false)
end

function Base.show(io::IO, mime::MIME"text/plain", obj::NCStatePolynomial)
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

function Base.string(obj::NCStatePolynomial)
    return join(
        map(zip(obj.coeffs, obj.nc_state_words)) do (coef, ncsw)
            string(coef) * " * " * string(ncsw)
        end,
        " + ",
    )
end

function print_object(io::IO, obj::NCStatePolynomial; multiline::Bool)
    return multiline ? print(io, string(obj)) : Base.show_default(io, obj)
end

function Base.:(==)(a::NCStatePolynomial, b::NCStatePolynomial)
    a.nc_state_words != b.nc_state_words && return false
    a.coeffs != b.coeffs && return false
    return true
end

function Base.hash(ncsp::NCStatePolynomial, u::UInt)
    return hash(hash.(ncsp.coeffs, u), hash.(ncsp.nc_state_words, u))
end

function Base.:(+)(a::NCStatePolynomial{T1}, b::NCStatePolynomial{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    return ncstatepoly(T[a.coeffs; b.coeffs], [a.nc_state_words; b.nc_state_words])
end

function Base.:(+)(a::NCStatePolynomial{T}, b::NCStateWord) where {T}
    return ncstatepoly(T[a.coeffs; one(T)], [a.nc_state_words; b])
end

function Base.:(-)(a::NCStatePolynomial{T1}, b::NCStatePolynomial{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    return ncstatepoly(T[a.coeffs; -one(T); b.coeffs], [a.nc_state_words; b.nc_state_words])
end

function Base.:(-)(a::NCStatePolynomial{T}, b::NCStateWord) where {T}
    return ncstatepoly(T[a.coeffs; -one(T)], [a.nc_state_words; b])
end

Base.one(::NCStatePolynomial{T,ST}) where {T,ST} = NCStatePolynomial([one(T)], [one(NCStateWord{ST})])
function Base.one(::Type{NCStatePolynomial{T}}) where {T}
    return NCStatePolynomial([one(T)], [one(NCStateWord)])
end

function Base.zero(::NCStatePolynomial{T}) where {T}
    return NCStatePolynomial([zero(T)], [one(NCStateWord)])
end

function Base.zero(::Type{NCStatePolynomial{T}}) where {T}
    return NCStatePolynomial([zero(T)], [one(NCStateWord)])
end

function variables(ncsp::NCStatePolynomial)
    return sorted_union(variables.(ncsp.nc_state_words)...)
end

function degree(ncsp::NCStatePolynomial)
    return reduce(max, degree.(ncsp.nc_state_words))
end

monomials(ncsp::NCStatePolynomial) = ncsp.nc_state_words

coefficients(ncsp::NCStatePolynomial) = ncsp.coeffs

terms(ncsp::NCStatePolynomial) = zip(ncsp.coeffs, ncsp.nc_state_words)
