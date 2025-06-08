get_dim(cons::VectorConstraint) = cons.set isa MOI.PositiveSemidefiniteConeSquare ? JuMP.shape(cons).side_dimension : JuMP.shape(cons).dims[1]

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
