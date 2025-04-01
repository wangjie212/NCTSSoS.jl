struct StateWord{V,M}
    state_monos::Vector{Monomial{V,M}}
    function StateWord(monos::Vector{Monomial{V,M}}) where {V,M}
        filter!(!isone, sort!(monos))
        new{V,M}(isempty(monos) ? [one(Monomial{V,M})] : monos)
    end
end

StateWord(vars::Vector{Variable{V,M}}) where {V,M} = StateWord(monomial.(vars))

DynamicPolynomials.variables(sw::StateWord) = union(variables.(sw.state_monos)...)
Base.show(io::IO, sw::StateWord) = print(io, join(map(x -> "<$(x)>", sw.state_monos), " * "))
Base.:(==)(a::StateWord, b::StateWord) = length(a.state_monos) == length(b.state_monos) && all(a.state_monos .== b.state_monos) # NOTE: need to guarantee it is always sorted
Base.hash(a::StateWord) = hash(a.state_monos)
Base.isless(a::StateWord, b::StateWord) = isless(a.state_monos, b.state_monos)
Base.:(*)(a::StateWord, b::StateWord) = StateWord([a.state_monos; b.state_monos])
Base.:(*)(coef::T, a::StateWord{V,M}) where {V,M,T} = StatePolynomial([coef], [a])
Base.one(a::StateWord) = StateWord([one(a.state_monos[1])])
DynamicPolynomials.degree(sw::StateWord) = mapreduce(degree, +, sw.state_monos; init=zero(Int))

struct NCStateWord{V,M}
    sw::StateWord{V,M}
    nc_word::Monomial{V,M}
end

Base.adjoint(a::NCStateWord{V,M}) where {V,M} = NCStateWord{V,M}(a.sw, star(a.nc_word))
Base.:(*)(a::NCStateWord{V,M}, b::NCStateWord{V,M}) where {V,M} = NCStateWord{V,M}(a.sw * b.sw, a.nc_word * b.nc_word)
Base.:(*)(coef::T, a::NCStateWord{V,M}) where {V,M,T} = StatePolynomialOp([coef * a.sw], [a.nc_word])
Base.:(==)(a::NCStateWord{V,M}, b::NCStateWord{V,M}) where {V,M} = a.sw == b.sw && a.nc_word == b.nc_word
Base.hash(a::NCStateWord) = hash((hash(a.sw), hash(a.nc_word)))

Base.show(io::IO, ncsw::NCStateWord) = print(io, string(ncsw.sw) * " ⋅ " * string(ncsw.nc_word))
function Base.isless(a::NCStateWord{V,M}, b::NCStateWord{V,M}) where {V,M}
    comp_val = compare(a.nc_word, b.nc_word)
    return comp_val < 0 || (comp_val == 0 && isless(a.sw, b.sw))
end
function expval(a::NCStateWord)
    StateWord([a.sw.state_monos; a.nc_word])
end

# FIXME: perhaps a better abstraction is NCStateWord -> StatePolynomialOp
# for the sake of progress, let's move on for now

struct StatePolynomial{V,M,T}
    coeffs::Vector{T}
    state_words::Vector{StateWord{V,M}}
    function StatePolynomial(coeffs::Vector{T}, state_words::Vector{StateWord{V,M}}) where {V,M,T}
        @assert length(coeffs) == length(state_words)
        uniq_state_words = sorted_unique(state_words)
        new{V,M,T}(map(x -> reduce(+, getindex.(Ref(coeffs), findall(==(x), state_words))), uniq_state_words), uniq_state_words)
    end
end

StatePolynomial(coeffs::Vector{T}, state_monos::Vector{Monomial{V,M}}) where {V,M,T} = StatePolynomial(coeffs, [StateWord([sm]) for sm in state_monos])

DynamicPolynomials.variables(sp::StatePolynomial) = union(variables.(sp.state_words)...)
DynamicPolynomials.degree(sp::StatePolynomial) = mapreduce(degree, max, sp.state_words)
DynamicPolynomials.monomials(sp::StatePolynomial) = sp.state_words
DynamicPolynomials.terms(sp::StatePolynomial) = collect(zip(sp.coeffs, sp.state_words))

Base.show(io::IO, sp::StatePolynomial) = print(io, join(map(x -> "$(x[1]) * $(x[2])", zip(sp.coeffs, sp.state_words)), " + "))
Base.:(==)(a::StatePolynomial, b::StatePolynomial) = (length(a.coeffs) == length(b.coeffs)) && mapfoldl(x -> (isequal(x[1], x[3]) && (x[2] == x[4])), &, zip(a.coeffs, a.state_words, b.coeffs, b.state_words))
Base.:(*)(a::StatePolynomial{V,M,T}, b::StatePolynomial{V,M,T}) where {V,M,T} =
    StatePolynomial(mapfoldl((args...) -> push!.(args...), product(zip(a.coeffs, a.state_words), zip(b.coeffs, b.state_words)); init=(T[], StateWord{V,M}[])) do ((ca, wa), (cb, wb))
        ca * cb, wa * wb
    end...)
Base.:(*)(n, a::StatePolynomial{V,M,T}) where {V,M,T} = StatePolynomial(T(n) .* a.coeffs, a.state_words)
Base.:(+)(a::StatePolynomial{V,M,T}, b::StatePolynomial{V,M,T}) where {V,M,T} = StatePolynomial([a.coeffs; b.coeffs], [a.state_words; b.state_words])
Base.one(a::StatePolynomial{V,M,T}) where {V,M,T} = StatePolynomial([one(T)], [one(a.state_words[1])])



# T: type of coefficient
# V: whether the variabels is NonCommutative{CreationOrder} or Commutative{CreationOrder}
# M: ordering of the monomials Graded{LexOrder} or else
# FIXME: inheriting from AbstractPolynomial gives StackOverflowError when trying to print in REPL
struct StatePolynomialOp{V,M,T}
    state_poly::Vector{StatePolynomial{V,M,T}}
    words::Vector{Monomial{V,M}}

    # TODO: get a more user friendly way of declaring StatePolynomial
    function StatePolynomialOp(state_polys::Vector{StatePolynomial{V,M,T}}, words::Vector{Monomial{V,M}}) where {V,M,T}
        @assert length(state_polys) == length(words) "Coefficients, formal words, and words must have the same length"
        uniq_words = sorted_unique(words)
        new{V,M,T}(map(x -> reduce(+, getindex.(Ref(state_polys), findall(==(x), words))), uniq_words), uniq_words)
    end
end

Base.show(io::IO, ncsp::StatePolynomialOp) = print(io, join(map(
        x -> (isone(x[2]) ? "$(string(x[1]))" : "$(string(x[1])) ⋅ $(string(x[2]))"),
        zip(ncsp.state_poly, ncsp.words)
    ), " + "))

Base.:(==)(a::StatePolynomialOp, b::StatePolynomialOp) = (a.state_poly == b.state_poly) && (a.words == b.words) # by constructor I alwasy guarantee no duplicate words and sorted
Base.:(+)(a::StatePolynomialOp, b::StatePolynomialOp) = StatePolynomialOp([a.state_poly; b.state_poly], [a.words; b.words])
Base.one(a::StatePolynomialOp) = StatePolynomialOp([one(a.state_poly[1])], [one(a.words[1])])
DynamicPolynomials.variables(ncsp::StatePolynomialOp) = sorted_union(variables.(ncsp.words)..., variables.(ncsp.state_poly)...)
DynamicPolynomials.degree(ncsp::StatePolynomialOp) = mapreduce(x -> sum(degree.(x)), max, zip(ncsp.state_poly, ncsp.words))
DynamicPolynomials.monomials(ncsp::StatePolynomialOp) =
    mapreduce(vcat, zip(ncsp.state_poly, ncsp.words)) do (sp, wd)
        [NCStateWord(c_word, wd) for c_word in monomials(sp)]
    end

function DynamicPolynomials.terms(ncsp::StatePolynomialOp)
    mapreduce(vcat, zip(ncsp.state_poly, ncsp.words)) do (sp, nc_word)
        [(sp_term[1], NCStateWord(sp_term[2], nc_word)) for sp_term in terms(sp)]
    end
end

function get_state_basis(variables::Vector{Variable{V,M}}, d::Int) where {V,M}
    return mapreduce(vcat, 0:d) do nc_deg
        nc_basis = monomials(variables, nc_deg)
        cw_deg = d - nc_deg
        cw_basis = unique!([
            begin
                interm = sort(filter(!isone, collect(c_word)))
                isempty(interm) ? [one(variables[1])] : interm
            end for c_word in product(ntuple(_ -> monomials(variables, 0:cw_deg), cw_deg)..., [one(variables[1])]) if sum(degree.(c_word)) <= cw_deg
        ])
        reshape(collect(product(cw_basis, nc_basis)), :)
    end
end
