function star(m::Monomial)
    return Monomial(reverse(m.vars), reverse(m.z))
end

function symmetric_canonicalize(monomial::Monomial)
    isempty(monomial.vars) && return monomial
    return min(monomial, star(monomial))
end

function symmetric_canonicalize(poly::Polynomial)
    return Polynomial(conj.(poly.coeffs), symmetric_canonicalize.(poly.monos))
end

# cyclic canonical is both cyclic and symmetric
function cyclic_canonicalize(monomial::Monomial)
    isempty(monomial.vars) && return monomial
    flatten_vars = mapreduce(
        idx -> fill(monomial.vars[idx], monomial.z[idx]), vcat, eachindex(monomial.z)
    )
    flatten_z = ones(Int, sum(monomial.z))
    return minimum(
        mapreduce(vcat, 1:sum(monomial.z)) do shift
            shifted_mono = Monomial(circshift!(flatten_vars, 1), circshift!(flatten_z, 1))
            [shifted_mono, star(shifted_mono)]
        end,
    )
end

function cyclic_canonicalize(poly::Polynomial)
    return Polynomial(poly.coeffs, cyclic_canonicalize.(poly.monos))
end

function get_basis(vars::Vector{Variable}, d::Int)
    return sort(
        mapreduce(vcat, 0:d) do dg
            vec(
                map(Iterators.product(repeat([vars], dg)...)) do cur_var
                    Monomial(cur_var, ones(dg))
                end,
            )
        end,
    )
end

function support(poly::Polynomial{T}, canonicalize::Function) where {T}
    return unique!(canonicalize.(poly.monos))
end

function neat_dot(x::Monomial, y::Monomial)
    return star(x) * y
end

# function neat_dot(x::NCStateWord{V,M}, y::NCStateWord{V,M}) where {V,M}
#     return adjoint(x) * y
# end

sorted_unique(xs) = sort(unique(xs))
sorted_union(xs...) = sort(union(xs...))

function _comm(mono::Monomial, comm_gps::Vector{Set{Variable}})
    map(comm_gps) do vars
        prod(zip(mono.vars, mono.z); init=Monomial([], [])) do (var, expon)
            var in vars ? var^expon : var^(zero(expon))
        end
    end
end

function _unipotent(mono::Monomial)
    isempty(mono.vars) && return mono
    prev_mono = mono
    local cur_mono
    while true
        cur_mono = prod(zip(prev_mono.vars, prev_mono.z); init=one(Monomial)) do (var, expo)
            var^(expo % 2)
        end
        cur_mono == prev_mono && break
        prev_mono = cur_mono
    end
    return cur_mono
end

# function _unipotent(ncsw::NCStateWord)
#     NCStateWord(_unipotent.(ncsw.sw), _unipotent(ncsw.nc_word))
# end

function _projective(mono::Monomial)
    prod(zip(mono.vars, mono.z); init=one(Monomial)) do (var, expo)
        var^(iszero(expo) ? expo : one(expo))
    end
end

# _projective(ncsw::NCStateWord) =
#     NCStateWord(_projective.(ncsw.sw), _projective(ncsw.nc_word))
