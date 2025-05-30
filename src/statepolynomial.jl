struct NCStateWord{V,M}
    sw::Vector{Monomial{V,M}}
    nc_word::Monomial{V,M}
    function NCStateWord(monos::Vector{Monomial{V,M}}, nc_word::Monomial{V,M}) where {V,M}
        filter!(!isone, sort!(monos))
        new{V,M}(isempty(monos) ? [one(Monomial{V,M})] : monos, nc_word)
    end
end

ς(m::Union{Monomial,Variable}) = NCStateWord([monomial(m)], one(m))

DynamicPolynomials.effective_variables(ncsw::NCStateWord) = union(effective_variables(ncsw.nc_word), effective_variables.(ncsw.sw)...)
DynamicPolynomials.degree(ncsw::NCStateWord) = degree(ncsw.nc_word) + sum(degree.(ncsw.sw))
# NOTE: if nc_word = one(x) where x is a dynamicpolynomial variable, it will appear in variables, this is potentially buggy
DynamicPolynomials.variables(ncsw::NCStateWord) = union(variables(ncsw.nc_word), variables.(ncsw.sw)...)

Base.adjoint(a::NCStateWord{V,M}) where {V,M} = NCStateWord(star.(a.sw), star(a.nc_word))
Base.:(*)(a::NCStateWord{V,M}, b::NCStateWord{V,M}) where {V,M} = NCStateWord([a.sw; b.sw], a.nc_word * b.nc_word)
Base.:(*)(a::NCStateWord{V,M}, b::Monomial{V,M}) where {V,M} = NCStateWord(a.sw, a.nc_word * b)
Base.:(*)(coef::T, a::NCStateWord{V,M}) where {V,M,T} = NCStateTerm(coef, a)
Base.:(==)(a::NCStateWord{V,M}, b::NCStateWord{V,M}) where {V,M} = (a.sw == b.sw) && a.nc_word == b.nc_word
Base.hash(a::NCStateWord) = hash((hash(a.sw), hash(a.nc_word)))
Base.one(::Type{NCStateWord{V,M}}) where {V,M} = NCStateWord([one(Monomial{V,M})], one(Monomial{V,M}))
Base.one(ncsw::NCStateWord) = NCStateWord([one(ncsw.sw[1])], one(ncsw.nc_word))
function Base.isless(a::NCStateWord{V,M}, b::NCStateWord{V,M}) where {V,M}
    comp_val = compare(a.nc_word, b.nc_word)
    return comp_val < 0 || (comp_val == 0 && isless(a.sw, b.sw))
end
Base.show(io::IO, ncsw::NCStateWord) = print(io, join(["<$(string(mono))>" for mono in ncsw.sw], " * ") * " * " * string(ncsw.nc_word))

expval(a::NCStateWord) = NCStateWord([a.sw; a.nc_word], one(a.nc_word))

struct NCStateTerm{V,M,T}
    coef::T
    ncstate_word::NCStateWord{V,M}
end

Base.:(==)(a::NCStateTerm, b::NCStateTerm) = isequal(a.coef, b.coef) && (a.ncstate_word == b.ncstate_word)
Base.:(+)(a::NCStateTerm, b::NCStateTerm) = NCStatePolynomial([a, b])
Base.:(-)(a::NCStateTerm{V,M,T}, b::NCStateTerm{V,M,T}) where {V,M,T} = NCStatePolynomial([a, -b])
Base.:(-)(a::NCStateTerm) = NCStateTerm(-a.coef, a.ncstate_word)
Base.:(*)(n, a::NCStateTerm{V,M,T}) where {V,M,T} = NCStateTerm(n * a.coef, a.ncstate_word)
Base.:(*)(a::NCStateTerm{V,M,T}, b::Monomial{V,M}) where {V,M,T} = NCStateTerm(a.coef, a.ncstate_word * b)
Base.:(*)(a::NCStateTerm{V,M,T}, b::NCStateTerm{V,M,T}) where {V,M,T} = NCStateTerm(a.coef * b.coef, a.ncstate_word * b.ncstate_word)
Base.hash(a::NCStateTerm) = hash((hash(a.coef), hash(a.ncstate_word)))
Base.one(::Type{NCStateTerm{V,M,T}}) where {V,M,T} = NCStateTerm(one(T), one(NCStateWord{V,M}))
Base.zero(::Type{NCStateTerm{V,M,T}}) where {V,M,T} = NCStateTerm(zero(T), one(NCStateWord{V,M}))
Base.isless(a::NCStateTerm, b::NCStateTerm) = isless(a.ncstate_word, b.ncstate_word)
DynamicPolynomials.degree(ncst::NCStateTerm) = degree(ncst.ncstate_word)
DynamicPolynomials.monomial(ncst::NCStateTerm) = ncst.ncstate_word
DynamicPolynomials.variables(ncst::NCStateTerm) = variables(ncst.ncstate_word)
function Base.show(io::IO, ncst::NCStateTerm)
    printstyled(io, ncst.coef; color=:yellow)
    print(io, " * " * string(ncst.ncstate_word))
end

# T: type of coefficient
# V: whether the variabels is NonCommutative{CreationOrder} or Commutative{CreationOrder}
# M: ordering of the monomials Graded{LexOrder} or else
struct NCStatePolynomial{V,M,T}
    nc_state_terms::Vector{NCStateTerm{V,M,T}}
    function NCStatePolynomial(nc_state_terms::Vector{NCStateTerm{V,M,T}}) where {V,M,T}
        uniq_words = sorted_unique([ncst.ncstate_word for ncst in nc_state_terms])
        new{V,M,T}(map(uniq_words) do word
            NCStateTerm(sum([ncst.coef for ncst in nc_state_terms[findall(x -> x.ncstate_word == (word), nc_state_terms)]]), word)
        end)
    end
end

Base.show(io::IO, ncsp::NCStatePolynomial) = print(io, join(string.(ncsp.nc_state_terms), " + "))

Base.:(==)(a::NCStatePolynomial, b::NCStatePolynomial) = all(a.nc_state_terms .== b.nc_state_terms) # by constructor I alwasy guarantee no duplicate words and sorted
Base.hash(ncsp::NCStatePolynomial) = hash(hash.(ncsp.nc_state_terms))
Base.:(*)(a::NCStatePolynomial, b::NCStatePolynomial) = NCStatePolynomial(vec(map(x -> x[1] * x[2], product(a.nc_state_terms, b.nc_state_terms))))
Base.:(*)(n, a::NCStatePolynomial{V,M,T}) where {V,M,T} = NCStatePolynomial(T(n) .* a.nc_state_terms)
Base.:(+)(a::NCStatePolynomial, b::NCStatePolynomial) = NCStatePolynomial([a.nc_state_terms; b.nc_state_terms])
Base.:(+)(a::NCStatePolynomial, b::NCStateTerm) = NCStatePolynomial([a.nc_state_terms; b])
Base.:(-)(a::NCStatePolynomial{V,M,T}, b::NCStatePolynomial{V,M,T}) where {V,M,T} = NCStatePolynomial([a.nc_state_terms; -one(T) .* b.nc_state_terms])
Base.:(-)(a::NCStatePolynomial{V,M,T}, b::NCStateTerm) where {V,M,T} = NCStatePolynomial([a.nc_state_terms; -one(T) * b])

Base.one(::NCStatePolynomial{V,M,T}) where {V,M,T} = NCStatePolynomial([one(NCStateTerm{V,M,T})])
Base.zero(::NCStatePolynomial{V,M,T}) where {V,M,T} = NCStatePolynomial([zero(NCStateTerm{V,M,T})])
Base.zero(::Type{NCStatePolynomial{V,M,T}}) where {V,M,T} = NCStatePolynomial([zero(NCStateTerm{V,M,T})])
DynamicPolynomials.variables(ncsp::NCStatePolynomial) = sorted_union(variables.(ncsp.nc_state_terms)...)
DynamicPolynomials.degree(ncsp::NCStatePolynomial) = reduce(max, degree.(ncsp.nc_state_terms))
DynamicPolynomials.monomials(ncsp::NCStatePolynomial) = monomial.(ncsp.nc_state_terms)
DynamicPolynomials.terms(ncsp::NCStatePolynomial) = ncsp.nc_state_terms

function get_state_basis(variables::Vector{Variable{V,M}}, d::Int, reducer) where {V,M}
    return map(a -> NCStateWord(a[1], a[2]), mapreduce(vcat, 0:d) do nc_deg
        nc_basis = unique(filter(m -> degree(m) == nc_deg, reducer.(monomials(variables, nc_deg))))
        cw_deg = d - nc_deg
        cw_basis = unique!([
            begin
                interm = sort(filter(!isone, collect(c_word)))
                isempty(interm) ? [one(variables[1])] : interm

            end for c_word in product(ntuple(_ -> unique!(symmetric_canonicalize.(reducer.(monomials(variables, 0:cw_deg)))), cw_deg)..., [one(variables[1])]) if sum(degree.(c_word)) <= cw_deg
        ])
        reshape(collect(product(cw_basis, nc_basis)), :)
    end)
end

for symb in [:symmetric_canonicalize, :cyclic_canonicalize]
    take_adj = (symb == :symmetric_canonicalize ? :adjoint : :identity)
    eval(quote
        function $(symb)(ncsw::NCStateWord)
            NCStateWord($(symb).(ncsw.sw), $(symb)(ncsw.nc_word))
        end

        function $(symb)(ncst::NCStateTerm)
            NCStateTerm($(take_adj)(ncst.coef), $(symb)(ncst.ncstate_word))
        end

        function $(symb)(spo::NCStatePolynomial)
            NCStatePolynomial($(symb).(spo.nc_state_terms))
        end
    end)
end
