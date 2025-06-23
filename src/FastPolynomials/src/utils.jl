"""
    sorted_unique(xs)

Returns sorted unique elements from a collection.

# Arguments
- `xs`: Collection to process

# Returns
- Sorted vector of unique elements
"""
sorted_unique(xs) = sort!(unique(xs))
sorted_unique!(xs) = sort!(unique!(xs))

"""
    sorted_union(xs...)

Returns sorted union of multiple collections.

# Arguments
- `xs...`: Variable number of collections

# Returns
- Sorted vector containing union of all input collections
"""
sorted_union(xs...) = sort(union(xs...))
