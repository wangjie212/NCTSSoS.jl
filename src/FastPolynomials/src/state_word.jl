abstract type StateType end
struct Arbitrary <: StateType end
struct MaxEntangled <: StateType end
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
struct StateWord{ST<:StateType}
    state_monos::Vector{Monomial}
    function StateWord{ST}(monos::Vector{Monomial}) where {ST}
        filter!(!isone, sort!(monos))
        return new{ST}(isempty(monos) ? [one(Monomial)] : monos)
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
ς(m::Union{Monomial,Variable}) = StateWord{Arbitrary}([monomial(m)])

"""
    tr(m::Union{Monomial,Variable})

Creates a StateWord with MaxEntangled state type from a monomial or variable.
The name tr stands for trace, commonly used in quantum contexts.

# Arguments
- `m::Union{Monomial,Variable}`: Monomial or variable to convert

# Returns
- `StateWord{MaxEntangled}`: StateWord with MaxEntangled state type containing the input as a single state monomial

# Example
```jldoctest; setup=:(using NCTSSoS.FastPolynomials; using NCTSSoS.FastPolynomials: tr)
julia> @ncpolyvar x y z
(x, y, z)

julia> tr(x^2*y)
tr(x²y¹)
```
"""
tr(m::Union{Monomial,Variable}) = StateWord{MaxEntangled}([monomial(m)])

variables(sw::StateWord) = sorted_union(variables.(sw.state_monos)...)

star(sw::StateWord{ST}) where {ST} = StateWord{ST}(star.(sw.state_monos))

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

function Base.string(obj::StateWord{ST}) where {ST}
    lencloser = ST == MaxEntangled ? "tr(" : "<"
    rencloser = ST == MaxEntangled ? ")" : ">"
    return join(map(x -> lencloser * "$(string(x))" * rencloser, obj.state_monos), " * ")
end

function print_object(io::IO, obj::StateWord; multiline::Bool)
    return multiline ? print(io, string(obj)) : Base.show_default(io, obj)
end

function Base.cmp(a::StateWord, b::StateWord)
    degree(a) != degree(b) && return degree(a) < degree(b) ? -1 : 1
    return cmp(a.state_monos, b.state_monos)
end

Base.:(==)(a::StateWord, b::StateWord) = iszero(cmp(a, b))

Base.hash(a::StateWord, u::UInt) = hash(a.state_monos, u)

Base.isless(a::StateWord, b::StateWord) = cmp(a, b) < 0

function Base.:(*)(a::StateWord{ST}, b::StateWord{ST}) where {ST}
    return StateWord{ST}([a.state_monos; b.state_monos])
end
Base.:(*)(a::StateWord{ST}, b::Monomial) where {ST} = NCStateWord{ST}(a, b)

Base.one(::Type{StateWord{ST}}) where {ST} = StateWord{ST}([one(Monomial)])
Base.one(_::StateWord{ST}) where {ST} = one(StateWord{ST})
Base.zero(::Type{StateWord{ST}}) where {ST} = 0.0 * one(StateWord{ST})
Base.zero(::StateWord{ST}) where {ST} = 0.0 * one(StateWord{ST})

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
<x²y¹> * x¹z¹
```
"""
struct NCStateWord{ST}
    sw::StateWord{ST}
    nc_word::Monomial
end

function NCStateWord(::Type{ST}, sw::Vector, nc_word) where {ST}
    return NCStateWord(StateWord{ST}(monomial.(sw)), monomial(nc_word))
end
NCStateWord(sw::StateWord{ST}) where {ST} = NCStateWord(ST, sw, one(Monomial))

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

function Base.one(::Type{NCStateWord{ST}}) where {ST}
    return NCStateWord(one(StateWord{ST}), one(Monomial))
end
Base.one(_::NCStateWord{ST}) where {ST} = one(NCStateWord{ST})
Base.zero(::Type{NCStateWord{ST}}) where {ST} = 0.0 * one(NCStateWord{ST})
Base.zero(::NCStateWord{ST}) where {ST} = 0.0 * one(NCStateWord{ST})

"""
    expval(a::NCStateWord)

Computes the expectation value by combining state monomials with the non-commutative word.

# Arguments
- `a::NCStateWord`: The NCStateWord

# Returns
- `StateWord`: StateWord containing all state monomials plus the non-commutative word
"""
expval(a::NCStateWord{ST}) where {ST} = StateWord{ST}([a.sw.state_monos; a.nc_word])
