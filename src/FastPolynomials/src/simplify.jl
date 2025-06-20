struct SimplifyAlgorithm
    comm_gps::Vector{Vector{Variable}}
    is_unipotent::Bool
    is_projective::Bool
end

function simplify(m::Monomial, sa::SimplifyAlgorithm)
    cxs = _comm(m, sa.comm_gps)
    return prod(
        sa.is_unipotent ? _unipotent.(cxs) : (sa.is_projective ? _projective.(cxs) : cxs)
    )
end

function simplify(sw::StateWord, sa::SimplifyAlgorithm)
    return StateWord(simplify.(sw.state_monos, sa))
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
