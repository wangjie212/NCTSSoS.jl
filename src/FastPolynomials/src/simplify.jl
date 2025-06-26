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

function simplify(sw::StateWord{ST}, sa::SimplifyAlgorithm) where ST
    return StateWord{ST}(simplify.(sw.state_monos, Ref(sa)))
end

function simplify(ncsw::NCStateWord, sa::SimplifyAlgorithm)
    return NCStateWord(simplify(ncsw.sw, sa), simplify(ncsw.nc_word, sa))
end

"""
    cyclic_canonicalize(mono::Monomial)

Canonicalizes a monomial using both cyclic and symmetric operations.
Finds the minimum among all cyclic shifts and their adjoints.
Chclic canonical is both cyclic and symmetric

# Arguments
- `mono::Monomial`: The mono to canonicalize

# Returns
- `Monomial`: Canonicalized monomial (minimum across all cyclic shifts and their stars)
"""
function cyclic_canonicalize(mono::Monomial, sa::SimplifyAlgorithm)
    isempty(mono.vars) && return mono
    flatten_vars = mapreduce(
        idx -> fill(mono.vars[idx], mono.z[idx]), vcat, eachindex(mono.z)
    )
    flatten_z = ones(Int, sum(mono.z))
    return minimum(
        mapreduce(vcat, 1:sum(mono.z)) do shift
            shifted_mono = monomial(circshift!(flatten_vars, 1), circshift!(flatten_z, 1))
            [simplify(shifted_mono,sa), simplify(star(shifted_mono),sa)]
        end,
    )
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

function canonicalize(sw::StateWord{MaxEntangled}, sa::SimplifyAlgorithm)
    return StateWord{MaxEntangled}(cyclic_canonicalize.(sw.state_monos, Ref(sa)))
end

function canonicalize(sw::StateWord{Arbitrary}, sa::SimplifyAlgorithm)
    return StateWord{Arbitrary}(symmetric_canonicalize.(sw.state_monos, Ref(sa)))
end

function canonicalize(ncsw::NCStateWord, sa::SimplifyAlgorithm)
    return NCStateWord(
        canonicalize(ncsw.sw, sa), symmetric_canonicalize(ncsw.nc_word, sa)
    )
end

"""
    canonicalize(poly::Polynomial)

Canonicalizes a polynomial by applying symmetric canonicalization to all monomials.

It only makes sense of symmetric canonicalize a polynomial or non-commuting
varaibles because trace polynomial is taken care at StatePoly

# Arguments
- `poly::Polynomial`: The polynomial to canonicalize

# Returns
- `Polynomial`: Canonicalized polynomial with conjugated coefficients and canonicalized monomials
"""
function canonicalize(poly::Polynomial, sa::SimplifyAlgorithm)
    return Polynomial(conj.(poly.coeffs), symmetric_canonicalize.(poly.monos, Ref(sa)))
end

function _unipotent(sw::StateWord{ST}) where ST
    return StateWord{ST}(_unipotent.(sw.state_monos))
end

function _unipotent(ncsw::NCStateWord)
    return NCStateWord(_unipotent(ncsw.sw), _unipotent(ncsw.nc_word))
end

_projective(sw::StateWord{ST}) where ST = StateWord{ST}(_projective.(sw.state_monos))

function _projective(ncsw::NCStateWord)
    return NCStateWord(_projective(ncsw.sw), _projective(ncsw.nc_word))
end

"""
ς(w) = ς(w') stated in https://arxiv.org/abs/2301.12513, Section 2.1
"""
function get_state_basis(
    ::Type{ST}, variables::Vector{Variable}, d::Int, sa::SimplifyAlgorithm
) where ST
    return unique!(
        map(
            a -> NCStateWord(StateWord{ST}(a[1]), a[2]),
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
    ::Type{NCStatePolynomial{T,ST}}, variables::Vector{Variable}, d::Int, sa::SimplifyAlgorithm
) where {T,ST}
    return get_state_basis(ST, variables, d, sa)
end


function canonicalize(sp::StatePolynomial,sa::SimplifyAlgorithm)
    return StatePolynomial((sp.coeffs), canonicalize.(sp.state_words, Ref(sa)))
end

function canonicalize(ncsp::NCStatePolynomial, sa::SimplifyAlgorithm)
    return NCStatePolynomial(
        (ncsp.coeffs), canonicalize.(ncsp.nc_state_words, Ref(sa))
    )
end

function get_basis(
    ::Type{Polynomial{T}}, variables::Vector{Variable}, d::Int, ::SimplifyAlgorithm
) where {T}
    return get_basis(variables, d)
end
