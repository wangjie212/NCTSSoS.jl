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
```jldoctest; setup:=(using NCTSSoS.FastPolynomials: ς)
"""
ς(m::Union{Monomial,Variable}) = StateWord([Monomial(m)])

"""
    variables(sw::StateWord)

Extracts all variables appearing in a StateWord.

# Arguments
- `sw::StateWord`: The StateWord

# Returns
- Set of all variables in the StateWord's monomials
"""
variables(sw::StateWord) = union(variables.(sw.state_monos)...)

"""
    degree(sw::StateWord)

Computes the total degree of a StateWord (sum of degrees of all state monomials).

# Arguments
- `sw::StateWord`: The StateWord

# Returns
- `Int`: Total degree of the StateWord
"""
function degree(sw::StateWord)
    return mapreduce(degree, +, sw.state_monos; init=zero(Int))
end

# 2-argument show, used by Array show, print(obj) and repr(obj), keep it short
function Base.show(io::IO, obj::StateWord)
    return print_object(io, obj; multiline=false)
end

# the 3-argument show used by display(obj) on the REPL
function Base.show(io::IO, mime::MIME"text/plain", obj::StateWord)
    # you can add IO options if you want
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

"""
    Base.string(obj::StateWord)

Converts a StateWord to its string representation with angle brackets.

# Arguments
- `obj::StateWord`: The StateWord to convert

# Returns
- `String`: String representation with each monomial in angle brackets, joined by " * "
"""
function Base.string(obj::StateWord)
    return join(map(x -> "<$(string(x))>", obj.state_monos), " * ")
end

"""
    print_object(io, obj; multiline)

Prints a StateWord object to an IO stream.

# Arguments
- `io::IO`: The output stream
- `obj::StateWord`: The StateWord to print
- `multiline::Bool`: Whether to use string representation or default format

# Returns
- Nothing (prints to IO stream)
"""
function print_object(io::IO, obj::StateWord; multiline::Bool)
    return multiline ? print(io, string(obj)) : Base.show_default(io, obj)
end

"""
    Base.cmp(a::StateWord, b::StateWord)

Compares two StateWords using graded lexicographic ordering.
First compares by total degree, then by state monomials vectors.

# Arguments
- `a::StateWord`: First StateWord
- `b::StateWord`: Second StateWord

# Returns
- `Int`: -1 if a < b, 0 if a == b, 1 if a > b
"""
function Base.cmp(a::StateWord, b::StateWord)
    degree(a) != degree(b) && return degree(a) < degree(b) ? -1 : 1
    return cmp(a.state_monos, b.state_monos)
end

"""
    Base.:(==)(a::StateWord, b::StateWord)

Tests equality between two StateWords.

# Arguments
- `a::StateWord`: First StateWord
- `b::StateWord`: Second StateWord

# Returns
- `Bool`: True if StateWords are equal
"""
Base.:(==)(a::StateWord, b::StateWord) = iszero(cmp(a, b))

"""
    Base.hash(a::StateWord, u::UInt)

Computes hash value for a StateWord based on its state monomials.

# Arguments
- `a::StateWord`: The StateWord to hash
- `u::UInt`: Hash seed value

# Returns
- `UInt`: Hash value of the state monomials (requires sorted invariant)
"""
# NOTE: need to guarantee it is always sorted
Base.hash(a::StateWord, u::UInt) = hash(a.state_monos, u)

"""
    Base.isless(a::StateWord, b::StateWord)

Determines if one StateWord is less than another using graded lexicographic ordering.

# Arguments
- `a::StateWord`: First StateWord
- `b::StateWord`: Second StateWord

# Returns
- `Bool`: True if a is less than b
"""
Base.isless(a::StateWord, b::StateWord) = cmp(a, b) < 0

"""
    Base.:(*)(a::StateWord, b::StateWord)

Multiplies two StateWords by concatenating their state monomials.

# Arguments
- `a::StateWord`: First StateWord
- `b::StateWord`: Second StateWord

# Returns
- `StateWord`: Product with concatenated state monomials
"""
Base.:(*)(a::StateWord, b::StateWord) = StateWord([a.state_monos; b.state_monos])

"""
    Base.:(*)(a::StateWord, b::Monomial)

Multiplies a StateWord by a monomial, creating an NCStateWord.

# Arguments
- `a::StateWord`: StateWord (commutative part)
- `b::Monomial`: Monomial (non-commutative part)

# Returns
- `NCStateWord`: Non-commutative state word with the StateWord and monomial
"""
Base.:(*)(a::StateWord, b::Monomial) = NCStateWord(a, b)
# Base.:(*)(coef::T, a::StateWord) where {T} = StatePolynomial([coef], [a])

"""
    Base.one(a::StateWord)

Returns the multiplicative identity for StateWord type (instance method).

# Arguments
- `a::StateWord`: StateWord instance

# Returns
- `StateWord`: Identity StateWord with single identity monomial
"""
Base.one(a::StateWord) = StateWord([one(a.state_monos[1])])

"""
    Base.one(::Type{StateWord})

Returns the multiplicative identity for StateWord type (type method).

# Arguments
- `::Type{StateWord}`: StateWord type

# Returns
- `StateWord`: Identity StateWord with single identity monomial
"""
Base.one(::Type{StateWord}) = StateWord([one(Monomial)])

"""
    NCStateWord

A non-commutative state word combining a commutative StateWord with a non-commutative Monomial.
Represents mixed commutative-noncommutative polynomial expressions.

# Fields
- `sw::StateWord`: Commutative state word part
- `nc_word::Monomial`: Non-commutative monomial part
"""
struct NCStateWord
    sw::StateWord
    nc_word::Monomial
end

"""
    degree(ncsw::NCStateWord)

Computes the total degree of an NCStateWord (sum of both parts).

# Arguments
- `ncsw::NCStateWord`: The NCStateWord

# Returns
- `Int`: Total degree of state word and non-commutative word
"""
degree(ncsw::NCStateWord) = degree(ncsw.nc_word) + degree(ncsw.sw)

"""
    variables(ncsw::NCStateWord)

Extracts all variables appearing in an NCStateWord.

# Arguments
- `ncsw::NCStateWord`: The NCStateWord

# Returns
- Set of all variables from both the state word and non-commutative word
"""
function variables(ncsw::NCStateWord)
    return union(variables(ncsw.nc_word), variables(ncsw.sw))
end

"""
    Base.adjoint(a::NCStateWord)

Computes the adjoint of an NCStateWord by taking the star of the non-commutative part.

# Arguments
- `a::NCStateWord`: The NCStateWord

# Returns
- `NCStateWord`: Adjoint with same state word and starred non-commutative word
"""
Base.adjoint(a::NCStateWord) = NCStateWord(a.sw, star(a.nc_word))

"""
    Base.:(*)(a::NCStateWord, b::NCStateWord)

Multiplies two NCStateWords by multiplying corresponding parts.

# Arguments
- `a::NCStateWord`: First NCStateWord
- `b::NCStateWord`: Second NCStateWord

# Returns
- `NCStateWord`: Product with multiplied state words and non-commutative words
"""
function Base.:(*)(a::NCStateWord, b::NCStateWord)
    return NCStateWord(a.sw * b.sw, a.nc_word * b.nc_word)
end
# Base.:(*)(coef::T, a::NCStateWord) where {T} = NCStatePolynomial([coef], [a])
#

"""
    Base.cmp(a::NCStateWord, b::NCStateWord)

Compares two NCStateWords using graded ordering with non-commutative word priority.
First compares by total degree, then by non-commutative word, then by state word.

# Arguments
- `a::NCStateWord`: First NCStateWord
- `b::NCStateWord`: Second NCStateWord

# Returns
- `Int`: -1 if a < b, 0 if a == b, 1 if a > b
"""
function Base.cmp(a::NCStateWord, b::NCStateWord)
    degree(a) != degree(b) && return degree(a) < degree(b) ? -1 : 1
    nc_word_res = cmp(a.nc_word, b.nc_word)
    iszero(nc_word_res) && return cmp(a.sw, b.sw)
    return nc_word_res
end

"""
    Base.isless(a::NCStateWord, b::NCStateWord)

Determines if one NCStateWord is less than another.

# Arguments
- `a::NCStateWord`: First NCStateWord
- `b::NCStateWord`: Second NCStateWord

# Returns
- `Bool`: True if a is less than b
"""
Base.isless(a::NCStateWord, b::NCStateWord) = cmp(a, b) < 0

"""
    Base.:(==)(a::NCStateWord, b::NCStateWord)

Tests equality between two NCStateWords.

# Arguments
- `a::NCStateWord`: First NCStateWord
- `b::NCStateWord`: Second NCStateWord

# Returns
- `Bool`: True if NCStateWords are equal
"""
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
    return print_object(io, obj; multiline=false)
end

# the 3-argument show used by display(obj) on the REPL
function Base.show(io::IO, mime::MIME"text/plain", obj::NCStateWord)
    # you can add IO options if you want
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

"""
    Base.string(obj::NCStateWord)

Converts an NCStateWord to its string representation.

# Arguments
- `obj::NCStateWord`: The NCStateWord to convert

# Returns
- `String`: String representation with state word and non-commutative word separated by space
"""
Base.string(obj::NCStateWord) = string(obj.sw) * " " * string(obj.nc_word)

"""
    print_object(io, obj; multiline)

Prints an NCStateWord object to an IO stream.

# Arguments
- `io::IO`: The output stream
- `obj::NCStateWord`: The NCStateWord to print
- `multiline::Bool`: Whether to use string representation or default format

# Returns
- Nothing (prints to IO stream)
"""
function print_object(io::IO, obj::NCStateWord; multiline::Bool)
    return multiline ? print(io, string(obj)) : Base.show_default(io, obj)
end

"""
    Base.one(::Type{NCStateWord})

Returns the multiplicative identity for NCStateWord type.

# Arguments
- `::Type{NCStateWord}`: NCStateWord type

# Returns
- `NCStateWord`: Identity with identity state word and identity monomial
"""
Base.one(::Type{NCStateWord}) = NCStateWord(one(StateWord), one(Monomial))

"""
    expval(a::NCStateWord)

Computes the expectation value by combining state monomials with the non-commutative word.

# Arguments
- `a::NCStateWord`: The NCStateWord

# Returns
- `StateWord`: StateWord containing all state monomials plus the non-commutative word
"""
expval(a::NCStateWord) = StateWord([a.sw.state_monos; a.nc_word])

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
        return new{T}(uniq_coeffs, uniq_state_words)
    end
end

"""
    variables(sp::StatePolynomial)

Extracts all variables appearing in a StatePolynomial.

# Arguments
- `sp::StatePolynomial`: The StatePolynomial

# Returns
- Set of all variables in the polynomial's state words
"""
variables(sp::StatePolynomial) = union(variables.(sp.state_words)...)

"""
    degree(sp::StatePolynomial)

Computes the maximum degree among all state words in the polynomial.

# Arguments
- `sp::StatePolynomial`: The StatePolynomial

# Returns
- `Int`: Maximum degree of any state word in the polynomial
"""
degree(sp::StatePolynomial) = mapreduce(degree, max, sp.state_words)

function Base.show(io::IO, obj::StatePolynomial)
    return print_object(io, obj; multiline=false)
end

function Base.show(io::IO, mime::MIME"text/plain", obj::StatePolynomial)
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

"""
    Base.string(obj::StatePolynomial)

Converts a StatePolynomial to its string representation.

# Arguments
- `obj::StatePolynomial`: The StatePolynomial to convert

# Returns
- `String`: String representation with terms joined by " + "
"""
function Base.string(obj::StatePolynomial)
    return join(
        map(zip(obj.coeffs, obj.state_words)) do (coef, sw)
            string(coef) * " * " * string(sw)
        end,
        " + ",
    )
end

"""
    print_object(io, obj; multiline)

Prints a StatePolynomial object to an IO stream.

# Arguments
- `io::IO`: The output stream
- `obj::StatePolynomial`: The StatePolynomial to print
- `multiline::Bool`: Whether to use string representation or default format

# Returns
- Nothing (prints to IO stream)
"""
function print_object(io::IO, obj::StatePolynomial; multiline::Bool)
    return multiline ? print(io, string(obj)) : Base.show_default(io, obj)
end

"""
    Base.:(==)(a::StatePolynomial, b::StatePolynomial)

Tests equality between two StatePolynomials.

# Arguments
- `a::StatePolynomial`: First StatePolynomial
- `b::StatePolynomial`: Second StatePolynomial

# Returns
- `Bool`: True if both state words and coefficients are identical
"""
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

"""
    Base.:(*)(a::StatePolynomial{T}, b::StatePolynomial{T}) where {T}

Multiplies two StatePolynomials using distributive property.

# Arguments
- `a::StatePolynomial{T}`: First StatePolynomial
- `b::StatePolynomial{T}`: Second StatePolynomial

# Returns
- `StatePolynomial{T}`: Product with all pairwise coefficient and state word products
"""
function Base.:(*)(a::StatePolynomial{T}, b::StatePolynomial{T}) where {T}
    return StatePolynomial(
        vec([ac * bc for (ac, bc) in Iterators.product(a.coeffs, b.coeffs)]),
        vec([asw * bsw for (asw, bsw) in Iterators.product(a.state_words, b.state_words)]),
    )
end

"""
    Base.:(*)(a::StatePolynomial{T}, b::Monomial) where {T}

Multiplies a StatePolynomial by a monomial, creating an NCStatePolynomial.

# Arguments
- `a::StatePolynomial{T}`: StatePolynomial
- `b::Monomial`: Monomial to multiply by

# Returns
- `NCStatePolynomial{T}`: Non-commutative state polynomial with each state word multiplied by the monomial
"""
function Base.:(*)(a::StatePolynomial{T}, b::Monomial) where {T}
    return NCStatePolynomial(a.coeffs, [sw * b for sw in a.state_words])
end

"""
    Base.:(*)(n, a::StatePolynomial{T}) where {T}

Multiplies a StatePolynomial by a scalar.

# Arguments
- `n`: Scalar multiplier
- `a::StatePolynomial{T}`: StatePolynomial

# Returns
- `StatePolynomial{T}`: StatePolynomial with all coefficients scaled by n
"""
function Base.:(*)(n, a::StatePolynomial{T}) where {T}
    return StatePolynomial(T(n) .* a.coeffs, a.state_words)
end

"""
    Base.:(+)(a::StatePolynomial{T1}, b::StatePolynomial{T2}) where {T1,T2}

Adds two StatePolynomials with potentially different coefficient types.

# Arguments
- `a::StatePolynomial{T1}`: First StatePolynomial
- `b::StatePolynomial{T2}`: Second StatePolynomial

# Returns
- `StatePolynomial`: Sum with promoted coefficient type
"""
function Base.:(+)(a::StatePolynomial{T1}, b::StatePolynomial{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    return StatePolynomial(T[a.coeffs; b.coeffs], [a.state_words; b.state_words])
end

"""
    Base.:(+)(a::StatePolynomial{T}, b::StateWord) where {T}

Adds a StateWord to a StatePolynomial.

# Arguments
- `a::StatePolynomial{T}`: StatePolynomial
- `b::StateWord`: StateWord to add with coefficient 1

# Returns
- `StatePolynomial{T}`: Sum with the StateWord added as a new term
"""
function Base.:(+)(a::StatePolynomial{T}, b::StateWord) where {T}
    return StatePolynomial([a.coeffs; one(T)], [a.state_words; b])
end

"""
    Base.one(::StatePolynomial{T}) where {T}

Returns the multiplicative identity for StatePolynomial.

# Arguments
- `::StatePolynomial{T}`: StatePolynomial instance

# Returns
- `StatePolynomial{T}`: Identity polynomial with coefficient 1 and identity state word
"""
Base.one(::StatePolynomial{T}) where {T} = StatePolynomial([one(T)], [one(StateWord)])

"""
    Base.zero(::StatePolynomial{T}) where {T}

Returns the additive identity for StatePolynomial.

# Arguments
- `::StatePolynomial{T}`: StatePolynomial instance

# Returns
- `StatePolynomial{T}`: Zero polynomial with coefficient 0 and identity state word
"""
function Base.zero(::StatePolynomial{T}) where {T}
    return StatePolynomial([zero(T)], [one(StateWord)])
end

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
        return new{T}(uniq_coeffs, uniq_nc_state_words)
    end
end

function Base.show(io::IO, obj::NCStatePolynomial)
    return print_object(io, obj; multiline=false)
end

# the 3-argument show used by display(obj) on the REPL
function Base.show(io::IO, mime::MIME"text/plain", obj::NCStatePolynomial)
    # you can add IO options if you want
    multiline = get(io, :multiline, true)
    return print_object(io, obj; multiline=multiline)
end

"""
    Base.string(obj::NCStatePolynomial)

Converts an NCStatePolynomial to its string representation.

# Arguments
- `obj::NCStatePolynomial`: The NCStatePolynomial to convert

# Returns
- `String`: String representation with terms joined by " + "
"""
function Base.string(obj::NCStatePolynomial)
    return join(
        map(zip(obj.coeffs, obj.nc_state_words)) do (coef, ncsw)
            string(coef) * " * " * string(ncsw)
        end,
        " + ",
    )
end

"""
    print_object(io, obj; multiline)

Prints an NCStatePolynomial object to an IO stream.

# Arguments
- `io::IO`: The output stream
- `obj::NCStatePolynomial`: The NCStatePolynomial to print
- `multiline::Bool`: Whether to use string representation or default format

# Returns
- Nothing (prints to IO stream)
"""
function print_object(io::IO, obj::NCStatePolynomial; multiline::Bool)
    return multiline ? print(io, string(obj)) : Base.show_default(io, obj)
end

"""
    Base.:(==)(a::NCStatePolynomial, b::NCStatePolynomial)

Tests equality between two NCStatePolynomials.

# Arguments
- `a::NCStatePolynomial`: First NCStatePolynomial
- `b::NCStatePolynomial`: Second NCStatePolynomial

# Returns
- `Bool`: True if both NC state words and coefficients are identical
"""
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

"""
    Base.:(+)(a::NCStatePolynomial{T1}, b::NCStatePolynomial{T2}) where {T1,T2}

Adds two NCStatePolynomials with potentially different coefficient types.

# Arguments
- `a::NCStatePolynomial{T1}`: First NCStatePolynomial
- `b::NCStatePolynomial{T2}`: Second NCStatePolynomial

# Returns
- `NCStatePolynomial`: Sum with promoted coefficient type
"""
function Base.:(+)(a::NCStatePolynomial{T1}, b::NCStatePolynomial{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    return NCStatePolynomial(T[a.coeffs; b.coeffs], [a.nc_state_words; b.nc_state_words])
end

"""
    Base.:(+)(a::NCStatePolynomial{T}, b::NCStateWord) where {T}

Adds an NCStateWord to an NCStatePolynomial.

# Arguments
- `a::NCStatePolynomial{T}`: NCStatePolynomial
- `b::NCStateWord`: NCStateWord to add with coefficient 1

# Returns
- `NCStatePolynomial{T}`: Sum with the NCStateWord added as a new term
"""
function Base.:(+)(a::NCStatePolynomial{T}, b::NCStateWord) where {T}
    return NCStatePolynomial(T[a.coeffs; one(T)], [a.nc_state_words; b])
end

"""
    Base.:(-)(a::NCStatePolynomial{T1}, b::NCStatePolynomial{T2}) where {T1,T2}

Subtracts two NCStatePolynomials with potentially different coefficient types.

# Arguments
- `a::NCStatePolynomial{T1}`: First NCStatePolynomial (minuend)
- `b::NCStatePolynomial{T2}`: Second NCStatePolynomial (subtrahend)

# Returns
- `NCStatePolynomial`: Difference with promoted coefficient type
"""
function Base.:(-)(a::NCStatePolynomial{T1}, b::NCStatePolynomial{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    return NCStatePolynomial(
        T[a.coeffs; -one(T); b.coeffs], [a.nc_state_words; b.nc_state_words]
    )
end

"""
    Base.:(-)(a::NCStatePolynomial{T}, b::NCStateWord) where {T}

Subtracts an NCStateWord from an NCStatePolynomial.

# Arguments
- `a::NCStatePolynomial{T}`: NCStatePolynomial (minuend)
- `b::NCStateWord`: NCStateWord (subtrahend)

# Returns
- `NCStatePolynomial{T}`: Difference with the NCStateWord subtracted
"""
function Base.:(-)(a::NCStatePolynomial{T}, b::NCStateWord) where {T}
    return NCStatePolynomial(T[a.coeffs; -one(T)], [a.nc_state_words; b])
end

"""
    Base.one(::NCStatePolynomial{T}) where {T}

Returns the multiplicative identity for NCStatePolynomial.

# Arguments
- `::NCStatePolynomial{T}`: NCStatePolynomial instance

# Returns
- `NCStatePolynomial{T}`: Identity polynomial with coefficient 1 and identity NC state word
"""
function Base.one(::NCStatePolynomial{T}) where {T}
    return NCStatePolynomial([one(T)], [one(NCStateWord)])
end

"""
    Base.zero(::NCStatePolynomial{T}) where {T}

Returns the additive identity for NCStatePolynomial (instance method).

# Arguments
- `::NCStatePolynomial{T}`: NCStatePolynomial instance

# Returns
- `NCStatePolynomial{T}`: Zero polynomial with coefficient 0 and identity NC state word
"""
function Base.zero(::NCStatePolynomial{T}) where {T}
    return NCStatePolynomial([zero(T)], [one(NCStateWord)])
end

"""
    Base.zero(::Type{NCStatePolynomial{T}}) where {T}

Returns the additive identity for NCStatePolynomial (type method).

# Arguments
- `::Type{NCStatePolynomial{T}}`: NCStatePolynomial type

# Returns
- `NCStatePolynomial{T}`: Zero polynomial with coefficient 0 and identity NC state word
"""
function Base.zero(::Type{NCStatePolynomial{T}}) where {T}
    return NCStatePolynomial([zero(T)], [one(NCStateWord)])
end

"""
    variables(ncsp::NCStatePolynomial)

Extracts all variables appearing in an NCStatePolynomial.

# Arguments
- `ncsp::NCStatePolynomial`: The NCStatePolynomial

# Returns
- Sorted union of all variables in the polynomial's NC state words
"""
function variables(ncsp::NCStatePolynomial)
    return sorted_union(variables.(ncsp.nc_state_words)...)
end

function degree(ncsp::NCStatePolynomial)
    return reduce(max, degree.(ncsp.nc_state_words))
end

function get_state_basis(variables::Vector{Variable}, d::Int, reducer)
    return map(
        a -> NCStateWord(StateWord(a[1]), a[2]),
        mapreduce(vcat, 0:d) do nc_deg
            nc_basis = reducer.(monomials(variables, nc_deg))
            cw_deg = d - nc_deg
            cw_basis = unique!([
                begin
                    interm = sort(filter(!isone, collect(c_word)))
                    isempty(interm) ? [one(variables[1])] : interm

                    # if it is cyclic_canonicalize.(reducer.) the state basis matches why ?
                end for c_word in Iterators.product(
                    ntuple(
                        _ -> unique!(
                            symmetric_canonicalize.(reducer.(get_basis(variables, cw_deg))),
                        ),
                        cw_deg,
                    )...,
                    [one(variables[1])],
                ) if sum(degree.(c_word)) <= cw_deg
            ])
            reshape(collect(Iterators.product(cw_basis, nc_basis)), :)
        end,
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
                return StatePolynomial($(symb).(sp.state_words))
            end

            function $(symb)(spo::NCStatePolynomial)
                return NCStatePolynomial($(symb).(spo.nc_state_words))
            end
        end,
    )
end
