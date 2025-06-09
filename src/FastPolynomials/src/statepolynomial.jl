struct StateWord
    state_monos::Vector{Monomial}
    function StateWord(monos::Vector{Monomial})
        filter!(!isone, sort!(monos))
        return new(isempty(monos) ? [one(Monomial)] : monos)
    end
end

ς(m::Union{Monomial,Variable}) = StateWord([Monomial(m)])

variables(sw::StateWord) = union(variables.(sw.state_monos)...)

function degree(sw::StateWord)
    return mapreduce(degree, +, sw.state_monos; init=zero(Int))
end

function Base.show(io::IO, obj::StateWord)
    return print_object(io, obj; multiline=false)
end

# the 3-argument show used by display(obj) on the REPL
function Base.show(io::IO, mime::MIME"text/plain", obj::StateWord)
    # you can add IO options if you want
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
# NOTE: need to guarantee it is always sorted
Base.hash(a::StateWord, u::UInt) = hash(a.state_monos, u)
Base.isless(a::StateWord, b::StateWord) = cmp(a, b) < 0

Base.:(*)(a::StateWord, b::StateWord) = StateWord([a.state_monos; b.state_monos])
Base.:(*)(a::StateWord, b::Monomial) = NCStateWord(a, b)
# Base.:(*)(coef::T, a::StateWord) where {T} = StatePolynomial([coef], [a])
Base.one(a::StateWord) = StateWord([one(a.state_monos[1])])
Base.one(::Type{StateWord}) = StateWord([one(Monomial)])

struct NCStateWord
    sw::StateWord
    nc_word::Monomial
end

degree(ncsw::NCStateWord) = degree(ncsw.nc_word) + degree(ncsw.sw)
function variables(ncsw::NCStateWord)
    return union(variables(ncsw.nc_word), variables(ncsw.sw))
end

Base.adjoint(a::NCStateWord) = NCStateWord(a.sw, star(a.nc_word))
function Base.:(*)(a::NCStateWord, b::NCStateWord)
    return NCStateWord(a.sw * b.sw, a.nc_word * b.nc_word)
end
# Base.:(*)(coef::T, a::NCStateWord) where {T} = NCStatePolynomial([coef], [a])
#

function Base.cmp(a::NCStateWord, b::NCStateWord)
    degree(a) != degree(b) && return degree(a) < degree(b) ? -1 : 1
    nc_word_res = cmp(a.nc_word, b.nc_word)
    iszero(nc_word_res) && return cmp(a.sw, b.sw)
    return nc_word_res
end

Base.isless(a::NCStateWord, b::NCStateWord) = cmp(a, b) < 0
#
Base.:(==)(a::NCStateWord, b::NCStateWord) = iszero(cmp(a, b))

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

Base.string(obj::NCStateWord) = string(obj.sw) * " " * string(obj.nc_word)

function print_object(io::IO, obj::NCStateWord; multiline::Bool)
    return multiline ? print(io, string(obj)) : Base.show_default(io, obj)
end

Base.one(::Type{NCStateWord}) = NCStateWord(one(StateWord), one(Monomial))

expval(a::NCStateWord) = StateWord([a.sw.state_monos; a.nc_word])

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

variables(sp::StatePolynomial) = union(variables.(sp.state_words)...)
degree(sp::StatePolynomial) = mapreduce(degree, max, sp.state_words)

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

function Base.:(*)(a::StatePolynomial{T}, b::Monomial) where {T}
    return NCStatePolynomial(a.coeffs, [sw * b for sw in a.state_words])
end

function Base.:(*)(n, a::StatePolynomial{T}) where {T}
    return StatePolynomial(T(n) .* a.coeffs, a.state_words)
end

function Base.:(+)(a::StatePolynomial{T1}, b::StatePolynomial{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    return StatePolynomial(T[a.coeffs; b.coeffs], [a.state_words; b.state_words])
end

function Base.:(+)(a::StatePolynomial{T}, b::StateWord) where {T}
    return StatePolynomial([a.coeffs; one(T)], [a.state_words; b])
end

Base.one(::StatePolynomial{T}) where {T} = StatePolynomial([one(T)], [one(StateWord)])
function Base.zero(::StatePolynomial{T}) where {T}
    return StatePolynomial([zero(T)], [one(StateWord)])
end

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

function Base.one(::NCStatePolynomial{T}) where {T}
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
