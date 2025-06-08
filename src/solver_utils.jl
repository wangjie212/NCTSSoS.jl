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
    flatten_vars = mapreduce(idx -> fill(monomial.vars[idx], monomial.z[idx]), vcat, eachindex(monomial.z))
    flatten_z = ones(Int, sum(monomial.z))
    return minimum(mapreduce(vcat, 1:sum(monomial.z)) do shift
        shifted_mono = Monomial(circshift!(flatten_vars, 1), circshift!(flatten_z, 1))
        [shifted_mono, star(shifted_mono)]
    end)
end

function cyclic_canonicalize(poly::Polynomial)
    return Polynomial(poly.coeffs, cyclic_canonicalize.(poly.monos))
end

function get_basis(vars::Vector{Variable}, d::Int)
    return sort(mapreduce(vcat, 0:d) do dg
        vec(map(Iterators.product(repeat([vars], dg)...)) do cur_var
            Monomial(cur_var, ones(dg))
        end)
    end)
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

get_dim(cons::VectorConstraint) = cons.set isa MOI.PositiveSemidefiniteConeSquare ? JuMP.shape(cons).side_dimension : JuMP.shape(cons).dims[1]

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

_projective(mono::Monomial) =
    prod(zip(mono.vars, mono.z);init=one(Monomial)) do (var, expo)
        var^(iszero(expo) ? expo : one(expo))
    end

# _projective(ncsw::NCStateWord) =
#     NCStateWord(_projective.(ncsw.sw), _projective(ncsw.nc_word))

function reducer(pop::PolyOpt)
    function (x)
        cxs = _comm(x, pop.comm_gps)
        return pop.is_unipotent ? _unipotent.(cxs) : (pop.is_projective ? _projective.(cxs) : cxs)
    end
end

# function reducer(spop::StatePolyOpt)
#     function (x)
#         cxs = _comm(x, spop.comm_gps)
#         return spop.is_unipotent ? _unipotent.(cxs) : (spop.is_projective ? _projective.(cxs) : cxs)
#     end
# end



# function _comm(ncsw::NCStateWord{V,M}, comm_gps::Vector{Set{Variable{V,M}}}) where {V,M}
#     [NCStateWord(prod.(_comm.(ncsw.sw, Ref(comm_gps))), prod(_comm(ncsw.nc_word, comm_gps)))]
# end


# function get_mom_matrix(mom_problem::MomentProblem)
#     _, mom_loc = findmax(get_dim, constraint_object.(mom_problem.constraints))
#     return value.(mom_problem.constraints[mom_loc])
# end

# function binary_search(a::T, b::Vector{T}) where {T}
#     left, right = 1, length(b)
#     while left <= right
#         mid = div(left + right, 2)
#         if b[mid] < a
#             left = mid + 1
#         elseif b[mid] > a
#             right = mid - 1
#         else
#             return mid
#         end
#     end
#     return 0
# end

# function get_mom_matrix(sos_problem::SOSProblem)
# end

# function minimizer_extraction(M::Matrix)

# end

# first do _comm and separate into different monomials
# then simplify with projective unipotent
# then do canonicalization
