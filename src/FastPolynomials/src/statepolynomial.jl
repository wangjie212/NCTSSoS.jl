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

function Base.:(==)(a::StateWord, b::StateWord)
    return length(a.state_monos) == length(b.state_monos) &&
           all(a.state_monos .== b.state_monos)
end
# NOTE: need to guarantee it is always sorted
Base.hash(a::StateWord) = hash(a.state_monos)
Base.isless(a::StateWord, b::StateWord) = isless(a.state_monos, b.state_monos)

Base.:(*)(a::StateWord, b::StateWord) = StateWord([a.state_monos; b.state_monos])
# Base.:(*)(a::StateWord, b::Monomial) = NCStateWord(a, b)
Base.:(*)(coef::T, a::StateWord) where {T} = StatePolynomial([coef], [a])
Base.one(a::StateWord) = StateWord([one(a.state_monos[1])])
Base.one(::Type{StateWord}) = StateWord([one(Monomial)])

# struct NCStateWord{V,M}
#     sw::StateWord{V,M}
#     nc_word::Monomial{V,M}
# end

# function DynamicPolynomials.effective_variables(ncsw::NCStateWord)
#     return union(effective_variables(ncsw.nc_word), effective_variables(ncsw.sw))
# end
# DynamicPolynomials.degree(ncsw::NCStateWord) = degree(ncsw.nc_word) + degree(ncsw.sw)
# function DynamicPolynomials.variables(ncsw::NCStateWord)
#     return union(variables(ncsw.nc_word), variables(ncsw.sw))
# end

# Base.adjoint(a::NCStateWord{V,M}) where {V,M} = NCStateWord{V,M}(a.sw, star(a.nc_word))
# function Base.:(*)(a::NCStateWord{V,M}, b::NCStateWord{V,M}) where {V,M}
#     return NCStateWord{V,M}(a.sw * b.sw, a.nc_word * b.nc_word)
# end
# Base.:(*)(coef::T, a::NCStateWord{V,M}) where {V,M,T} = NCStateTerm(coef, a)
# function Base.:(==)(a::NCStateWord{V,M}, b::NCStateWord{V,M}) where {V,M}
#     return a.sw == b.sw && a.nc_word == b.nc_word
# end
# Base.hash(a::NCStateWord) = hash((hash(a.sw), hash(a.nc_word)))
# function Base.show(io::IO, ncsw::NCStateWord)
#     return print(io, string(ncsw.sw) * " ⋅ " * string(ncsw.nc_word))
# end
# function Base.one(::Type{NCStateWord{V,M}}) where {V,M}
#     return NCStateWord{V,M}(one(StateWord{V,M}), one(Monomial{V,M}))
# end
# function Base.isless(a::NCStateWord{V,M}, b::NCStateWord{V,M}) where {V,M}
#     comp_val = compare(a.nc_word, b.nc_word)
#     return comp_val < 0 || (comp_val == 0 && isless(a.sw, b.sw))
# end

# expval(a::NCStateWord) = StateWord([a.sw.state_monos; a.nc_word])

struct StatePolynomial{T}
    coeffs::Vector{T}
    state_words::Vector{StateWord}
    function StatePolynomial(coeffs::Vector{T}, state_words::Vector{StateWord}) where {T}
        uniq_state_words = sorted_unique([st.state_word for st in state_terms])
        return new{T}(
            map(uniq_state_words) do sw
                StateTerm(
                    sum([
                        st.coef for st in
                        getindex.(
                            Ref(state_terms), findall(x -> x.state_word == sw, state_terms)
                        )
                    ]),
                    sw,
                )
            end,
        )
    end
end

DynamicPolynomials.variables(sp::StatePolynomial) = union(variables.(sp.state_terms)...)
DynamicPolynomials.degree(sp::StatePolynomial) = mapreduce(degree, max, sp.state_terms)
DynamicPolynomials.monomials(sp::StatePolynomial) = [st.state_word for st in sp.state_terms]
DynamicPolynomials.terms(sp::StatePolynomial) = sp.state_terms

Base.show(io::IO, sp::StatePolynomial) = print(io, join(string.(sp.state_terms), " + "))
Base.:(==)(a::StatePolynomial, b::StatePolynomial) = all(a.state_terms .== b.state_terms)
Base.hash(a::StatePolynomial) = hash(hash.(a.state_terms))
function Base.:(*)(a::StatePolynomial{V,M,T}, b::StatePolynomial{V,M,T}) where {V,M,T}
    return StatePolynomial(
        vec(map(x -> x[1] * x[2], product(a.state_terms, b.state_terms)))
    )
end
function Base.:(*)(a::StatePolynomial{V,M,T}, b::Monomial{V,M}) where {V,M,T}
    return StatePolynomialOp([st * b for st in a.state_terms])
end
function Base.:(*)(n, a::StatePolynomial{V,M,T}) where {V,M,T}
    return StatePolynomial(T(n) .* a.state_terms)
end
function Base.:(+)(a::StatePolynomial{V,M,T}, b::StatePolynomial{V,M,T}) where {V,M,T}
    return StatePolynomial([a.state_terms; b.state_terms])
end
function Base.:(+)(a::StatePolynomial{V,M,T}, b::StateTerm{V,M,T}) where {V,M,T}
    return StatePolynomial([a.state_terms; b])
end
Base.one(::StatePolynomial{V,M,T}) where {V,M,T} = StatePolynomial([one(StateTerm{V,M,T})])
function Base.zero(::StatePolynomial{V,M,T}) where {V,M,T}
    return StatePolynomial([zero(StateTerm{V,M,T})])
end

# # T: type of coefficient
# # V: whether the variabels is NonCommutative{CreationOrder} or Commutative{CreationOrder}
# # M: ordering of the monomials Graded{LexOrder} or else
# # FIXME: inheriting from AbstractPolynomial gives StackOverflowError when trying to print in REPL
# struct StatePolynomialOp{V,M,T}
#     nc_state_terms::Vector{NCStateTerm{V,M,T}}
#     # TODO: get a more user friendly way of declaring StatePolynomial
#     function StatePolynomialOp(nc_state_terms::Vector{NCStateTerm{V,M,T}}) where {V,M,T}
#         uniq_words = sorted_unique([ncst.ncstate_word for ncst in nc_state_terms])
#         return new{V,M,T}(
#             map(uniq_words) do word
#                 NCStateTerm(
#                     sum([
#                         ncst.coef for ncst in nc_state_terms[findall(
#                             x -> x.ncstate_word == (word), nc_state_terms
#                         )]
#                     ]),
#                     word,
#                 )
#             end,
#         )
#     end
# end

# function Base.show(io::IO, ncsp::StatePolynomialOp)
#     return print(io, join(string.(ncsp.nc_state_terms), " + "))
# end

# function Base.:(==)(a::StatePolynomialOp, b::StatePolynomialOp)
#     return all(a.nc_state_terms .== b.nc_state_terms)
# end # by constructor I alwasy guarantee no duplicate words and sorted
# Base.hash(ncsp::StatePolynomialOp) = hash(hash.(ncsp.nc_state_terms))
# function Base.:(+)(a::StatePolynomialOp, b::StatePolynomialOp)
#     return StatePolynomialOp([a.nc_state_terms; b.nc_state_terms])
# end
# Base.:(+)(a::StatePolynomialOp, b::NCStateTerm) = StatePolynomialOp([a.nc_state_terms; b])
# function Base.:(-)(a::StatePolynomialOp{V,M,T}, b::StatePolynomialOp{V,M,T}) where {V,M,T}
#     return StatePolynomialOp([a.nc_state_terms; -one(T) .* b.nc_state_terms])
# end
# function Base.:(-)(a::StatePolynomialOp{V,M,T}, b::NCStateTerm) where {V,M,T}
#     return StatePolynomialOp([a.nc_state_terms; -one(T) * b])
# end

# function Base.one(::StatePolynomialOp{V,M,T}) where {V,M,T}
#     return StatePolynomialOp([one(NCStateTerm{V,M,T})])
# end
# function Base.zero(::StatePolynomialOp{V,M,T}) where {V,M,T}
#     return StatePolynomialOp([zero(NCStateTerm{V,M,T})])
# end
# function Base.zero(::Type{StatePolynomialOp{V,M,T}}) where {V,M,T}
#     return StatePolynomialOp([zero(NCStateTerm{V,M,T})])
# end
# function DynamicPolynomials.variables(ncsp::StatePolynomialOp)
#     return sorted_union(variables.(ncsp.nc_state_terms)...)
# end
# function DynamicPolynomials.degree(ncsp::StatePolynomialOp)
#     return reduce(max, degree.(ncsp.nc_state_terms))
# end
# DynamicPolynomials.monomials(ncsp::StatePolynomialOp) = monomial.(ncsp.nc_state_terms)
# DynamicPolynomials.terms(ncsp::StatePolynomialOp) = ncsp.nc_state_terms

# function get_state_basis(variables::Vector{Variable{V,M}}, d::Int, reducer) where {V,M}
#     return map(
#         a -> NCStateWord(StateWord(a[1]), a[2]),
#         mapreduce(vcat, 0:d) do nc_deg
#             nc_basis = reducer.(monomials(variables, nc_deg))
#             cw_deg = d - nc_deg
#             cw_basis = unique!([
#                 begin
#                     interm = sort(filter(!isone, collect(c_word)))
#                     isempty(interm) ? [one(variables[1])] : interm

#                     # if it is cyclic_canonicalize.(reducer.) the state basis matches why ?
#                 end for c_word in product(
#                     ntuple(
#                         _ -> unique!(
#                             symmetric_canonicalize.(
#                                 reducer.(monomials(variables, 0:cw_deg))
#                             ),
#                         ),
#                         cw_deg,
#                     )...,
#                     [one(variables[1])],
#                 ) if sum(degree.(c_word)) <= cw_deg
#             ])
#             reshape(collect(product(cw_basis, nc_basis)), :)
#         end,
#     )
# end

# for symb in [:symmetric_canonicalize, :cyclic_canonicalize]
#     take_adj = (symb == :symmetric_canonicalize ? :adjoint : :identity)
#     eval(
#         quote
#             function $(symb)(sw::StateWord)
#                 return StateWord($(symb).(sw.state_monos))
#             end

#             function $(symb)(ncsw::NCStateWord)
#                 return NCStateWord($(symb)(ncsw.sw), $(symb)(ncsw.nc_word))
#             end

#             function $(symb)(sp::StatePolynomial)
#                 return StatePolynomial($(symb).(sp.state_terms))
#             end

#             function $(symb)(spo::StatePolynomialOp)
#                 return StatePolynomialOp($(symb).(spo.nc_state_terms))
#             end
#         end,
#     )
# end
