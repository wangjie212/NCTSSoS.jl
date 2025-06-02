Base.cmp(a::Variable, b::Variable) = a.name == b.name ? 0 : (a.name < b.name ? -1 : 1)

Base.isless(a::Variable, b::Variable) = cmp(a, b) < 0

Base.in(a::Variable, collection::Vector{Variable}) =
	searchsortedfirst(collection, a, lt=cmp) != 0

function Base.cmp(a::Monomial, b::Monomial)
	degree(a) != degree(b) && return degree(a) < degree(b) ? -1 : 1
    first_nonzero_diff_idx = findfirst(
        i -> (a.vars[i] != b.vars[i]) && (!iszero(a.z[i] || !iszero(b.z[i]))),
        1:length(a.vars),
    )
    a_tuple, b_tuple = map(
        m -> ntuple(
            x -> m.z[x],
            isnothing(first_nonzero_diff_idx) ? length(a.z) : first_nonzero_diff_idx,
        ),
        (a, b),
    )
	return cmp(a_tuple, b_tuple)
end
