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

Base.one(::Type{StateWord}) = StateWord([one(Monomial)])
Base.one(_::StateWord) = one(StateWord)

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
