"""
    StateWord

A state word representing a product of state monomials.
Automatically removes identity elements and maintains sorted order.
A state monomial is the expectation value of a monomials of noncommutative
operators with respect to a state.

# Fields
- `state_monos::Vector{Monomial}`: Sorted vector of non-identity state monomials

# Constructor
Creates a StateWord from monomials, filtering out identity elements and sorting.
If all monomials are identity, stores a single identity monomial.
"""
struct StateWord
    state_monos::Vector{Monomial}
    function StateWord(monos::Vector{Monomial})
        filter!(!isone, sort!(monos))
        return new(isempty(monos) ? [one(Monomial)] : monos)
    end
end

"""
    ς(m::Union{Monomial,Variable})

Creates a StateWord from a monomial or variable using the Greek letter ς (sigma).

# Arguments
- `m::Union{Monomial,Variable}`: Monomial or variable to convert

# Returns
- `StateWord`: StateWord containing the input as a single state monomial

# Example
```jldoctest; setup=:(using NCTSSoS.FastPolynomials; using NCTSSoS.FastPolynomials: ς)
julia> @ncpolyvar x y z
(x, y, z)

julia> ς(x^2*y)
<x²y¹>
```
"""
ς(m::Union{Monomial,Variable}) = StateWord([Monomial(m)])

variables(sw::StateWord) = sorted_union(variables.(sw.state_monos)...)

star(sw::StateWord) = StateWord(star.(sw.state_monos))

function degree(sw::StateWord)
    return mapreduce(degree, +, sw.state_monos; init=zero(Int))
end

function Base.show(io::IO, obj::StateWord)
    return print_object(io, obj; multiline=true)
end

function Base.show(io::IO, mime::MIME"text/plain", obj::StateWord)
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

function Base.string(obj::StateWord)
    return join(map(x -> "<$(string(x))>", obj.state_monos), " * ")
end

function print_object(io::IO, obj::StateWord; multiline::Bool)
    return multiline ? print(io, string(obj)) : Base.show_default(io, obj)
end

function Base.cmp(a::StateWord, b::StateWord)
    degree(a) != degree(b) && return degree(a) < degree(b) ? -1 : 1
    return cmp(a.state_monos, b.state_monos)
end

Base.:(==)(a::StateWord, b::StateWord) = iszero(cmp(a, b))

"""
    Base.hash(a::StateWord, u::UInt)

Computes hash value for a StateWord based on its state monomials.
Need to guarantee it is always sorted

# Arguments
- `a::StateWord`: The StateWord to hash
- `u::UInt`: Hash seed value

# Returns
- `UInt`: Hash value of the state monomials (requires sorted invariant)
"""
Base.hash(a::StateWord, u::UInt) = hash(a.state_monos, u)

Base.isless(a::StateWord, b::StateWord) = cmp(a, b) < 0

Base.:(*)(a::StateWord, b::StateWord) = StateWord([a.state_monos; b.state_monos])

Base.:(*)(a::StateWord, b::Monomial) = NCStateWord(a, b)
Base.:(*)(coef::T, a::StateWord) where {T<:Number} = StatePolynomial([coef], [a])

Base.one(a::StateWord) = StateWord([one(a.state_monos[1])])
Base.one(::Type{StateWord}) = StateWord([one(Monomial)])

"""
    NCStateWord

A non-commutative state word combining a commutative StateWord with a non-commutative Monomial.
Represents mixed commutative-noncommutative polynomial expressions.

# Fields
- `sw::StateWord`: Commutative state word part
- `nc_word::Monomial`: Non-commutative monomial part

# Examples
```jldoctest; setup=:(using NCTSSoS.FastPolynomials; using NCTSSoS.FastPolynomials:NCStateWord, ς)
julia> @ncpolyvar x y z
(x, y, z)

julia> sw = ς(x^2*y)
<x²y¹>

julia> ncsw = sw * (x*z)
<x²y¹> x¹z¹
```
"""
struct NCStateWord
    sw::StateWord
    nc_word::Monomial
end

NCStateWord(sw::Vector, nc_word) = NCStateWord(StateWord(Monomial.(sw)), nc_word)
NCStateWord(sw::StateWord) = NCStateWord(sw, one(Monomial))

degree(ncsw::NCStateWord) = degree(ncsw.nc_word) + degree(ncsw.sw)

function variables(ncsw::NCStateWord)
    return sorted_union(variables(ncsw.nc_word), variables(ncsw.sw))
end

Base.adjoint(a::NCStateWord) = NCStateWord(star(a.sw), star(a.nc_word))

function neat_dot(x::NCStateWord, y::NCStateWord)
    return adjoint(x) * y
end

function Base.:(*)(a::NCStateWord, b::NCStateWord)
    return NCStateWord(a.sw * b.sw, a.nc_word * b.nc_word)
end

function Base.:(*)(a::StateWord, b::NCStateWord)
    return NCStateWord(a) * b
end

Base.:(*)(coef::Number, a::NCStateWord) = NCStatePolynomial([coef], [a])

function Base.cmp(a::NCStateWord, b::NCStateWord)
    degree(a) != degree(b) && return degree(a) < degree(b) ? -1 : 1
    nc_word_res = cmp(a.nc_word, b.nc_word)
    iszero(nc_word_res) && return cmp(a.sw, b.sw)
    return nc_word_res
end

Base.isless(a::NCStateWord, b::NCStateWord) = cmp(a, b) < 0

Base.:(==)(a::NCStateWord, b::NCStateWord) = iszero(cmp(a, b))

"""
    Base.hash(a::NCStateWord, u::UInt)

Computes hash value for an NCStateWord based on both components.

# Arguments
- `a::NCStateWord`: The NCStateWord to hash
- `u::UInt`: Hash seed value

# Returns
- `UInt`: Hash value combining hashes of state word and non-commutative word
"""
Base.hash(a::NCStateWord, u::UInt) = hash((hash(a.sw, u), hash(a.nc_word, u)))

function Base.show(io::IO, obj::NCStateWord)
    return print_object(io, obj; multiline=true)
end

function Base.show(io::IO, mime::MIME"text/plain", obj::NCStateWord)
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

Base.string(obj::NCStateWord) = string(obj.sw) * " * " * string(obj.nc_word)

function print_object(io::IO, obj::NCStateWord; multiline::Bool)
    return multiline ? print(io, string(obj)) : Base.show_default(io, obj)
end

Base.one(::Type{NCStateWord}) = NCStateWord(one(StateWord), one(Monomial))
Base.one(_::NCStateWord) = one(NCStateWord)

"""
    expval(a::NCStateWord)

Computes the expectation value by combining state monomials with the non-commutative word.

# Arguments
- `a::NCStateWord`: The NCStateWord

# Returns
- `StateWord`: StateWord containing all state monomials plus the non-commutative word
"""
expval(a::NCStateWord) = StateWord([a.sw.state_monos; a.nc_word])

function _unipotent(sw::StateWord)
    return StateWord(_unipotent.(sw.state_monos))
end

function _unipotent(ncsw::NCStateWord)
    return NCStateWord(_unipotent(ncsw.sw), _unipotent(ncsw.nc_word))
end

_projective(sw::StateWord) = StateWord(_projective.(sw.state_monos))

function _projective(ncsw::NCStateWord)
    return NCStateWord(_projective(ncsw.sw), _projective(ncsw.nc_word))
end

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
struct StatePolynomial{T}
    coeffs::Vector{T}
    state_words::Vector{StateWord}
    function StatePolynomial(coeffs::Vector{T}, state_words::Vector{StateWord}) where {T}
        uniq_state_words = sorted_unique(state_words)
        uniq_coeffs = zeros(T, length(uniq_state_words))
        for (coef, sw) in zip(coeffs, state_words)
            idx = searchsortedfirst(uniq_state_words, sw)
            uniq_coeffs[idx] += coef
        end
        nz_idcs = filter(a -> !iszero(uniq_coeffs[a]), 1:length(uniq_coeffs))
        return new{T}(uniq_coeffs[nz_idcs], uniq_state_words[nz_idcs])
    end
end

variables(sp::StatePolynomial) = sorted_union(variables.(sp.state_words)...)

degree(sp::StatePolynomial) = mapreduce(degree, max, sp.state_words)

monomials(sp::StatePolynomial) = sp.state_words

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

"""
    Base.hash(a::StatePolynomial, u::UInt)

Computes hash value for a StatePolynomial.

# Arguments
- `a::StatePolynomial`: The StatePolynomial to hash
- `u::UInt`: Hash seed value

# Returns
- `UInt`: Hash value combining coefficients and state words
"""
Base.hash(a::StatePolynomial, u::UInt) = hash(hash.(a.coeffs, u), hash.(a.state_words, u))

function Base.:(*)(a::StatePolynomial{T}, b::StatePolynomial{T}) where {T}
    return StatePolynomial(
        vec([ac * bc for (ac, bc) in Iterators.product(a.coeffs, b.coeffs)]),
        vec([asw * bsw for (asw, bsw) in Iterators.product(a.state_words, b.state_words)]),
    )
end

function Base.:(*)(a::StatePolynomial{T}, b::Monomial) where {T}
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
struct NCStatePolynomial{T}
    coeffs::Vector{T}
    nc_state_words::Vector{NCStateWord}

    function NCStatePolynomial(
        coeffs::Vector{T}, nc_state_words::Vector{NCStateWord}
    ) where {T}
        uniq_nc_state_words = sorted_unique(nc_state_words)
        uniq_coeffs = zeros(T, length(uniq_nc_state_words))
        for (coef, sw) in zip(coeffs, nc_state_words)
            idx = searchsortedfirst(uniq_nc_state_words, sw)
            uniq_coeffs[idx] += coef
        end
        nz_idcs = filter(a -> !iszero(uniq_coeffs[a]), 1:length(uniq_coeffs))
        return new{T}(uniq_coeffs[nz_idcs], uniq_nc_state_words[nz_idcs])
    end
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

"""
    Base.hash(ncsp::NCStatePolynomial, u::UInt)

Computes hash value for an NCStatePolynomial.

# Arguments
- `ncsp::NCStatePolynomial`: The NCStatePolynomial to hash
- `u::UInt`: Hash seed value

# Returns
- `UInt`: Hash value combining coefficients and NC state words
"""
function Base.hash(ncsp::NCStatePolynomial, u::UInt)
    return hash(hash.(ncsp.coeffs, u), hash.(ncsp.nc_state_words, u))
end

function Base.:(+)(a::NCStatePolynomial{T1}, b::NCStatePolynomial{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    return NCStatePolynomial(T[a.coeffs; b.coeffs], [a.nc_state_words; b.nc_state_words])
end

function Base.:(+)(a::NCStatePolynomial{T}, b::NCStateWord) where {T}
    return NCStatePolynomial(T[a.coeffs; one(T)], [a.nc_state_words; b])
end

function Base.:(-)(a::NCStatePolynomial{T1}, b::NCStatePolynomial{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    return NCStatePolynomial(
        T[a.coeffs; -one(T); b.coeffs], [a.nc_state_words; b.nc_state_words]
    )
end

function Base.:(-)(a::NCStatePolynomial{T}, b::NCStateWord) where {T}
    return NCStatePolynomial(T[a.coeffs; -one(T)], [a.nc_state_words; b])
end

Base.one(::NCStatePolynomial{T}) where {T} = NCStatePolynomial([one(T)], [one(NCStateWord)])

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

terms(ncsp::NCStatePolynomial) = zip(ncsp.coeffs, ncsp.nc_state_words)

function get_state_basis(variables::Vector{Variable}, d::Int, reducer)
    return reducer.(
        map(
            a -> NCStateWord(StateWord(a[1]), a[2]),
            mapreduce(vcat, 0:d) do nc_deg
                nc_basis = monomials(variables, nc_deg)
                cw_deg = d - nc_deg
                cw_basis = unique!([
                    begin
                        interm = sort(filter(!isone, collect(c_word)))
                        isempty(interm) ? [one(variables[1])] : interm
                    end for c_word in Iterators.product(
                        ntuple(
                            _ -> unique!(
                                symmetric_canonicalize.(get_basis(variables, cw_deg))
                            ),
                            cw_deg,
                        )...,
                        [one(variables[1])],
                    ) if sum(degree.(c_word)) <= cw_deg
                ])
                reshape(collect(Iterators.product(cw_basis, nc_basis)), :)
            end,
        ),
    )
end

for symb in [:symmetric_canonicalize, :cyclic_canonicalize]
    take_adj = (symb == :symmetric_canonicalize ? :adjoint : :identity)
    eval(
        quote
            function $(symb)(sw::StateWord)
                return StateWord($(symb).(sw.state_monos))
            end

            function $(symb)(ncsw::NCStateWord)
                return NCStateWord($(symb)(ncsw.sw), $(symb)(ncsw.nc_word))
            end

            function $(symb)(sp::StatePolynomial)
                return StatePolynomial((sp.coeffs), $(symb).(sp.state_words))
            end

            function $(symb)(ncsp::NCStatePolynomial)
                return NCStatePolynomial((ncsp.coeffs), $(symb).(ncsp.nc_state_words))
            end
        end,
    )
end
