# NOTE: for ncpolyvar generating monomials from `monomials` may give you x^0*y^2*z^0*x^1 terms
# which is not equal to y^2*x^1
# plue this bug: https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/118#issue-1412618512
# I will keep this function here
function remove_zero_degree(m::Monomial)
    isconstant(m) && return m
    return prod([x[1]^x[2] for x in filter(!(iszero ∘ last), collect(zip(m.vars, m.z)))])
end

function star(m::Monomial)
    return prod([x[1]^x[2] for x in filter(!(iszero ∘ last), reverse(collect(zip(m.vars, m.z))))])
end

function symmetric_canonicalize(monomial::Monomial)
    return min(monomial, star(monomial))
end

function symmetric_canonicalize(poly::Polynomial)
    return mapreduce(p -> coefficient(p)' * symmetric_canonicalize(monomial(p)), +, terms(poly); init=zero(poly))
end

function cyclic_canonicalize(monomial::Monomial)
    isconstant(monomial) && return monomial
    return minimum(mapreduce(vcat, 0:length(monomial.vars)-1) do shift
        shifted_mono = prod([var^expo for (var, expo) in circshift!(collect(zip(monomial.vars, monomial.z)), shift)])
        [shifted_mono, star(shifted_mono)]
    end)
end

function cyclic_canonicalize(poly::Polynomial)
    return mapreduce(p -> coefficient(p) * cyclic_canonicalize(monomial(p)), +, terms(poly); init=zero(poly))
end

function get_basis(vars::Vector{Variable{V,M}}, d::Int) where {V,M}
    # need to remove zero degree other wise sortting fails
    # return mapreduce(cur_d -> remove_zero_degree.(monomials(vars, cur_d)), vcat, 0:d)
    return remove_zero_degree.(sort(monomials(vars, 0:d)))
end

function support(poly::Polynomial{V,M,T}, canonicalize::Function) where {V,M,T}
    return canonicalize.(monomials(poly))
end

function neat_dot(x::Monomial{V,M}, y::Monomial{V,M}) where {V,M}
    # NOTE: the `*` in DynamicPolynomials sometimes creates monomials with degree 0, which we don't want
    return remove_zero_degree(star(x) * y)
end

function neat_dot(x::NCStateWord{V,M}, y::NCStateWord{V,M}) where {V,M}
    return adjoint(x) * y
end

sorted_unique(xs) = sort(unique(xs))
sorted_union(xs...) = sort(union(xs...))

get_dim(cons::VectorConstraint) = cons.set isa MOI.PositiveSemidefiniteConeSquare ? JuMP.shape(cons).side_dimension : JuMP.shape(cons).dims[1]

function _comm(mono::Monomial{V,M}, comm_gp::Set{Variable{V,M}}) where {V,M}
    return [prod(zip(mono.vars, mono.z); init=one(mono)) do (var, expo)
            var in comm_gp ? var^expo : var^(zero(expo))
        end,
        prod(zip(mono.vars, mono.z); init=one(mono)) do (var, expo)
            !(var in comm_gp) ? var^expo : var^(zero(expo))
        end]
end

function _unipotent(mono::Monomial)
    isconstant(mono) && return mono
    prev_mono = mono
    local cur_mono
    while true
        cur_mono = prod(zip(prev_mono.vars, prev_mono.z)) do (var, expo)
            var^(expo % 2)
        end
        cur_mono == prev_mono && break
        prev_mono = cur_mono
    end
    return cur_mono
end

function _unipotent(ncsw::NCStateWord)
    NCStateWord(_unipotent.(ncsw.sw), _unipotent(ncsw.nc_word))
end

_projective(mono::Monomial) =
    prod(zip(mono.vars, mono.z)) do (var, expo)
        var^(iszero(expo) ? expo : one(expo))
    end
_projective(ncsw::NCStateWord) =
    NCStateWord(_projective.(ncsw.sw), _projective(ncsw.nc_word))

function reducer(pop::PolyOpt)
    function (x)
        cxs = _comm(x, pop.comm_gp)
        return pop.is_unipotent ? _unipotent.(cxs) : (pop.is_projective ? _projective.(cxs) : cxs)
    end
end

function reducer(spop::StatePolyOpt)
    function (x)
        cxs = _comm(x, spop.comm_gps)
        return spop.is_unipotent ? _unipotent.(cxs) : (spop.is_projective ? _projective.(cxs) : cxs)
    end
end

function _comm(mono::Monomial{V,M}, comm_gps::Vector{Set{Variable{V,M}}}) where {V,M}
    map(comm_gps) do vars
        prod(zip(mono.vars, mono.z); init=first(vars)^0) do (var, expon)
            var in vars ? var^expon : var^(zero(expon))
        end
    end
end


function _comm(ncsw::NCStateWord{V,M}, comm_gps::Vector{Set{Variable{V,M}}}) where {V,M}
    [NCStateWord(prod.(_comm.(ncsw.sw, Ref(comm_gps))), prod(_comm(ncsw.nc_word, comm_gps)))]
end


function get_mom_matrix(mom_problem::MomentProblem)
    _, mom_loc = findmax(get_dim, constraint_object.(mom_problem.constraints))
    return value.(mom_problem.constraints[mom_loc])
end

function binary_search(a::T, b::Vector{T}) where T
    left, right = 1, length(b)
    while left <= right
        mid = div(left + right, 2)
        if b[mid] < a
            left = mid + 1
        elseif b[mid] > a
            right = mid - 1
        else
            return mid
        end
    end
    return 0
end

function get_mom_matrix(sos_problem::SOSProblem)
end

function minimizer_extraction(M::Matrix)

end

# first do _comm and separate into different monomials
# then simplify with projective unipotent
# then do canonicalization