@kwdef struct SimplifyAlgorithm
    comm_gps::Vector{Vector{Variable}} = Vector{Variable}[]
    is_unipotent::Bool = false
    is_projective::Bool = false
end

function simplify(m::Monomial, sa::SimplifyAlgorithm)
    cxs = _comm(m, sa.comm_gps)
    return prod(
        sa.is_unipotent ? _unipotent.(cxs) : (sa.is_projective ? _projective.(cxs) : cxs)
    )
end

function simplify(sw::StateWord, sa::SimplifyAlgorithm)
    return StateWord(simplify.(sw.state_monos, Ref(sa)))
end

function simplify(ncsw::NCStateWord, sa::SimplifyAlgorithm)
    return NCStateWord(simplify(ncsw.sw, sa), simplify(ncsw.nc_word, sa))
end

"""
    symmetric_canonicalize(mono::Monomial)

Canonicalizes a mono by taking the minimum between itself and its adjoint.

# Arguments
- `mono::Monomial`: The mono to canonicalize

# Returns
- `Monomial`: Canonicalized monomial (minimum of original and its star)
"""
function symmetric_canonicalize(mono::Monomial, sa::SimplifyAlgorithm)
    isempty(mono.vars) && return mono
    return min(simplify(mono, sa), simplify(star(mono), sa))
end

function symmetric_canonicalize(sw::StateWord, sa::SimplifyAlgorithm)
    return StateWord(symmetric_canonicalize.(sw.state_monos, Ref(sa)))
end

function symmetric_canonicalize(ncsw::NCStateWord, sa::SimplifyAlgorithm)
    return NCStateWord(
        symmetric_canonicalize(ncsw.sw, sa), symmetric_canonicalize(ncsw.nc_word, sa)
    )
end

"""
    symmetric_canonicalize(poly::Polynomial)

Canonicalizes a polynomial by applying symmetric canonicalization to all monomials.

# Arguments
- `poly::Polynomial`: The polynomial to canonicalize

# Returns
- `Polynomial`: Canonicalized polynomial with conjugated coefficients and canonicalized monomials
"""
function symmetric_canonicalize(poly::Polynomial, sa::SimplifyAlgorithm)
    return Polynomial(conj.(poly.coeffs), symmetric_canonicalize.(poly.monos, Ref(sa)))
end

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
ς(w) = ς(w') stated in https://arxiv.org/abs/2301.12513, Section 2.1
"""
function get_state_basis(variables::Vector{Variable}, d::Int, sa::SimplifyAlgorithm)
    return unique!(
        map(
            a -> NCStateWord(StateWord(a[1]), a[2]),
            mapreduce(vcat, 0:d) do nc_deg
                nc_basis = simplify.(monomials(variables, Val(nc_deg)), Ref(sa))
                cw_deg = d - nc_deg
                cw_basis = unique!([
                    begin
                        interm = sort(filter(!isone, collect(c_word)))
                        isempty(interm) ? [one(variables[1])] : interm
                    end for c_word in Iterators.product(
                        ntuple(
                            _ -> unique!(simplify.(get_basis(variables, cw_deg), Ref(sa))),
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

function get_basis(
    ::Type{NCStatePolynomial{T}}, variables::Vector{Variable}, d::Int, sa::SimplifyAlgorithm
) where {T}
    return get_state_basis(variables, d, sa)
end

for symb in [:symmetric_canonicalize, :cyclic_canonicalize]
    eval(
        quote
            function $(symb)(sp::StatePolynomial, sa::SimplifyAlgorithm)
                return StatePolynomial((sp.coeffs), $(symb).(sp.state_words, Ref(sa)))
            end

            function $(symb)(ncsp::NCStatePolynomial, sa::SimplifyAlgorithm)
                return NCStatePolynomial(
                    (ncsp.coeffs), $(symb).(ncsp.nc_state_words, Ref(sa))
                )
            end
        end,
    )
end

function get_basis(
    ::Type{Polynomial{T}}, variables::Vector{Variable}, d::Int, ::SimplifyAlgorithm
) where {T}
    return get_basis(variables, d)
end
