Base.cmp(a::Variable, b::Variable) = a.name == b.name ? 0 : (a.name < b.name ? -1 : 1)

Base.isless(a::Variable, b::Variable) = cmp(a, b) < 0

Base.in(a::Variable, collection::Vector{Variable}) =
	searchsortedfirst(collection, a, lt=cmp) != 0