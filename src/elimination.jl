# implement NoElimination on clique_decomp
struct NoElimination <: EliminationAlgorithm end
struct AsIsElimination <: EliminationAlgorithm end

function cliquetree(graph, alg::NoElimination, snd::SupernodeType)
    return cliquetree(complete_graph(nv(graph)), BFS(), snd)
end

"""
    clique_decomp(G::SimpleGraph, clique_alg::EliminationAlgorithm)

Decomposes a graph into cliques using the specified elimination algorithm.

# Arguments
- `G::SimpleGraph`: Input graph to decompose
- `clique_alg::EliminationAlgorithm`: Algorithm for clique tree elimination

# Returns
- `Vector{Vector{Int}}`: Vector of cliques, each containing vertex indices
"""
function clique_decomp(G::SimpleGraph, clique_alg::EliminationAlgorithm)
    label, tree = cliquetree(G, alg=clique_alg)
    return map(x -> label[x], collect(Vector{Int}, tree))
end

# TODO: https://github.com/wangjie212/NCTSSoS.jl/issues/45
