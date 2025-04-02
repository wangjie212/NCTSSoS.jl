struct StateWord{V,M}
    state_monos::Vector{Monomial{V,M}}
    function StateWord(monos::Vector{Monomial{V,M}}) where {V,M}
        filter!(!isone, sort!(monos))
        new{V,M}(isempty(monos) ? [one(Monomial{V,M})] : monos)
    end
end

DynamicPolynomials.variables(sw::StateWord) = union(variables.(sw.state_monos)...)
DynamicPolynomials.degree(sw::StateWord) = mapreduce(degree, +, sw.state_monos; init=zero(Int))
Base.show(io::IO, sw::StateWord) = print(io, join(map(x -> "<$(x)>", sw.state_monos), " * "))
Base.:(==)(a::StateWord, b::StateWord) = length(a.state_monos) == length(b.state_monos) && all(a.state_monos .== b.state_monos) # NOTE: need to guarantee it is always sorted
Base.hash(a::StateWord) = hash(a.state_monos)
Base.isless(a::StateWord, b::StateWord) = isless(a.state_monos, b.state_monos)
Base.:(*)(a::StateWord, b::StateWord) = StateWord([a.state_monos; b.state_monos])
Base.:(*)(coef::T, a::StateWord{V,M}) where {V,M,T} = StateTerm(coef, a)
Base.one(a::StateWord) = StateWord([one(a.state_monos[1])])
Base.one(::Type{StateWord{V,M}}) where {V,M} = StateWord([one(Monomial{V,M})])

struct NCStateWord{V,M}
    sw::StateWord{V,M}
    nc_word::Monomial{V,M}
end

DynamicPolynomials.degree(ncsw::NCStateWord) = degree(ncsw.nc_word) + degree(ncsw.sw)
DynamicPolynomials.variables(ncsw::NCStateWord) = union(variables(ncsw.nc_word), variables(ncsw.sw))
Base.adjoint(a::NCStateWord{V,M}) where {V,M} = NCStateWord{V,M}(a.sw, star(a.nc_word))
Base.:(*)(a::NCStateWord{V,M}, b::NCStateWord{V,M}) where {V,M} = NCStateWord{V,M}(a.sw * b.sw, a.nc_word * b.nc_word)
Base.:(*)(coef::T, a::NCStateWord{V,M}) where {V,M,T} = NCStateTerm(coef, a)
Base.:(==)(a::NCStateWord{V,M}, b::NCStateWord{V,M}) where {V,M} = a.sw == b.sw && a.nc_word == b.nc_word
Base.hash(a::NCStateWord) = hash((hash(a.sw), hash(a.nc_word)))
Base.show(io::IO, ncsw::NCStateWord) = print(io, string(ncsw.sw) * " ⋅ " * string(ncsw.nc_word))
Base.one(::Type{NCStateWord{V,M}}) where {V,M} = NCStateWord{V,M}(one(StateWord{V,M}), one(Monomial{V,M}))
function Base.isless(a::NCStateWord{V,M}, b::NCStateWord{V,M}) where {V,M}
    comp_val = compare(a.nc_word, b.nc_word)
    return comp_val < 0 || (comp_val == 0 && isless(a.sw, b.sw))
end
function expval(a::NCStateWord)
    StateWord([a.sw.state_monos; a.nc_word])
end


struct StateTerm{V,M,T}
    coef::T
    state_word::StateWord{V,M}
    function StateTerm(coef::T, state_word::StateWord{V,M}) where {V,M,T}
        new{V,M,T}(coef, state_word)
    end
end

DynamicPolynomials.variables(st::StateTerm) = variables(st.state_word)
DynamicPolynomials.degree(st::StateTerm) = degree(st.state_word)
Base.show(io::IO, st::StateTerm) = print(io, string(st.coef) * " * " * string(st.state_word))
Base.:(==)(a::StateTerm, b::StateTerm) = isequal(a.coef, b.coef) && (a.state_word == b.state_word)
Base.hash(a::StateTerm) = hash((hash(a.coef), hash(a.state_word)))
Base.:(*)(a::StateTerm, b::StateTerm) = StateTerm(a.coef * b.coef, a.state_word * b.state_word)
Base.:(*)(n, a::StateTerm{V,M,T}) where {V,M,T} = StateTerm(T(n) * a.coef, a.state_word)
Base.:(+)(a::StateTerm, b::StateTerm) = StatePolynomial([a, b])
Base.one(::Type{StateTerm{V,M,T}}) where {V,M,T} = StateTerm(one(T), one(StateWord{V,M}))
Base.zero(::Type{StateTerm{V,M,T}}) where {V,M,T} = StateTerm(zero(T), one(StateWord{V,M}))

struct NCStateTerm{V,M,T}
    coef::T
    ncstate_word::NCStateWord{V,M}
    function NCStateTerm(coef::T, ncstate_word::NCStateWord{V,M}) where {V,M,T}
        new{V,M,T}(coef, ncstate_word)
    end
end

Base.:(==)(a::NCStateTerm, b::NCStateTerm) = isequal(a.coef, b.coef) && (a.ncstate_word == b.ncstate_word)
Base.:(+)(a::NCStateTerm, b::NCStateTerm) = StatePolynomialOp([a.coef, b.coef],[a.ncstate_word,b.ncstate_word])
Base.hash(a::NCStateTerm) = hash((hash(a.coef), hash(a.ncstate_word)))
Base.show(io::IO, ncst::NCStateTerm) = print(io, string(ncst.coef) * " * " * string(ncst.ncstate_word))
Base.one(::Type{NCStateTerm{V,M,T}}) where {V,M,T} = NCStateTerm(one(T), one(NCStateWord{V,M}))
Base.zero(::Type{NCStateTerm{V,M,T}}) where {V,M,T} = NCStateTerm(zero(T), one(NCStateWord{V,M}))
Base.isless(a::NCStateTerm, b::NCStateTerm) = isless(a.ncstate_word, b.ncstate_word)
DynamicPolynomials.degree(ncst::NCStateTerm) = degree(ncst.ncstate_word)
DynamicPolynomials.monomial(ncst::NCStateTerm) = ncst.ncstate_word
DynamicPolynomials.variables(ncst::NCStateTerm) = variables(ncst.ncstate_word)

struct StatePolynomial{V,M,T}
    state_terms::Vector{StateTerm{V,M,T}}
    function StatePolynomial(state_terms::Vector{StateTerm{V,M,T}}) where {V,M,T}
        uniq_state_words = sorted_unique([st.state_word for st in state_terms])
        new{V,M,T}(map(uniq_state_words) do sw
            StateTerm(sum([st.coef for st in getindex.(Ref(state_terms), findall(x -> x.state_word == sw, state_terms))]), sw)
        end)
    end
end

DynamicPolynomials.variables(sp::StatePolynomial) = union(variables.(sp.state_terms)...)
DynamicPolynomials.degree(sp::StatePolynomial) = mapreduce(degree, max, sp.state_terms)
DynamicPolynomials.monomials(sp::StatePolynomial) = [st.state_word for st in sp.state_terms]
DynamicPolynomials.terms(sp::StatePolynomial) = sp.state_terms

Base.show(io::IO, sp::StatePolynomial) = print(io, join(string.(sp.state_terms), " + "))
Base.:(==)(a::StatePolynomial, b::StatePolynomial) = all(a.state_terms .== b.state_terms)
Base.hash(a::StatePolynomial) = hash(hash.(a.state_terms))
Base.:(*)(a::StatePolynomial{V,M,T}, b::StatePolynomial{V,M,T}) where {V,M,T} = StatePolynomial(vec(map(x -> x[1] * x[2], product(a.state_terms, b.state_terms))))
Base.:(*)(n, a::StatePolynomial{V,M,T}) where {V,M,T} = StatePolynomial(T(n) .* a.state_terms)
Base.:(+)(a::StatePolynomial{V,M,T}, b::StatePolynomial{V,M,T}) where {V,M,T} = StatePolynomial([a.state_terms; b.state_terms])
Base.one(::StatePolynomial{V,M,T}) where {V,M,T} = StatePolynomial([one(StateTerm{V,M,T})])
Base.zero(::StatePolynomial{V,M,T}) where {V,M,T} = StatePolynomial([zero(StateTerm{V,M,T})])



# T: type of coefficient
# V: whether the variabels is NonCommutative{CreationOrder} or Commutative{CreationOrder}
# M: ordering of the monomials Graded{LexOrder} or else
# FIXME: inheriting from AbstractPolynomial gives StackOverflowError when trying to print in REPL
struct StatePolynomialOp{V,M,T}
    nc_state_terms::Vector{NCStateTerm{V,M,T}}
    # TODO: get a more user friendly way of declaring StatePolynomial
    function StatePolynomialOp(nc_state_terms::Vector{NCStateTerm{V,M,T}}) where {V,M,T}
        uniq_words = sorted_unique([ncst.ncstate_word for ncst in nc_state_terms])
        new{V,M,T}(map(uniq_words) do word
            NCStateTerm(sum([ncst.coef for ncst in nc_state_terms[findall(x -> x.ncstate_word == (word), nc_state_terms)]]), word)
        end)
    end
end

Base.show(io::IO, ncsp::StatePolynomialOp) = print(io, join(string.(ncsp.nc_state_terms), " + "))

Base.:(==)(a::StatePolynomialOp, b::StatePolynomialOp) = all(a.nc_state_terms .== b.nc_state_terms) # by constructor I alwasy guarantee no duplicate words and sorted
Base.hash(ncsp::StatePolynomialOp) = hash(hash.(ncsp.nc_state_terms))
Base.:(+)(a::StatePolynomialOp, b::StatePolynomialOp) = StatePolynomialOp([a.nc_state_terms; b.nc_state_terms])
Base.:(+)(a::StatePolynomialOp,b::NCStateTerm) = StatePolynomialOp([a.nc_state_terms; b])
Base.one(::StatePolynomialOp{V,M,T}) where {V,M,T} = StatePolynomialOp([one(NCStateTerm{V,M,T})])
Base.zero(::StatePolynomialOp{V,M,T}) where {V,M,T} = StatePolynomialOp([zero(NCStateTerm{V,M,T})])
Base.zero(::Type{StatePolynomialOp{V,M,T}}) where {V,M,T} = StatePolynomialOp([zero(NCStateTerm{V,M,T})])
DynamicPolynomials.variables(ncsp::StatePolynomialOp) = sorted_union(variables.(ncsp.nc_state_terms)...)
DynamicPolynomials.degree(ncsp::StatePolynomialOp) = reduce(max, degree.(ncsp.nc_state_terms))
DynamicPolynomials.monomials(ncsp::StatePolynomialOp) = monomial.(ncsp.nc_state_terms)
DynamicPolynomials.terms(ncsp::StatePolynomialOp) = ncsp.nc_state_terms

function get_state_basis(variables::Vector{Variable{V,M}}, d::Int) where {V,M}
    return map(a -> NCStateWord(StateWord(a[1]), a[2]), mapreduce(vcat, 0:d) do nc_deg
        nc_basis = monomials(variables, nc_deg)
        cw_deg = d - nc_deg
        cw_basis = unique!([
            begin
                interm = sort(filter(!isone, collect(c_word)))
                isempty(interm) ? [one(variables[1])] : interm
            end for c_word in product(ntuple(_ -> monomials(variables, 0:cw_deg), cw_deg)..., [one(variables[1])]) if sum(degree.(c_word)) <= cw_deg
        ])
        reshape(collect(product(cw_basis, nc_basis)), :)
    end)
end
