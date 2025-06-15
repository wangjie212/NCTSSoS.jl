"""
    CorrelativeSparsity

Structure representing the correlative sparsity pattern of a polynomial optimization problem.

# Fields
- `cliques::Vector{Vector{Variable}}`: Groups of variables that form cliques in the sparsity graph
- `cliques_cons::Vector{Vector{Int}}`: Constraint indices assigned to each clique
- `global_cons::Vector{Int}`: Constraint indices not captured by any single clique
- `cliques_idcs_bases::Vector{Vector{Vector{Monomial}}}`: Monomial bases for indexing moment/localizing matrices within each clique
"""
struct CorrelativeSparsity
    cliques::Vector{Vector{Variable}}
    cliques_cons::Vector{Vector{Int}}
    # FIXME: add test case for difference
    global_cons::Vector{Int}
    cliques_idcs_bases::Vector{Vector{Vector{Monomial}}}
end

function show(io::IO, cs::CorrelativeSparsity)
    max_size = maximum(length.(cs.cliques))
    println(io, "Correlative Sparsity: \n")
    println(io, "   maximum size: $max_size")
    for clique_i in 1:length(cs.cliques)
        println(io, "   Clique $clique_i: ")
        println(io, "       Variables: ", cs.cliques[clique_i])
        println(io, "       Constraints: ")
        for con_j in 1:length(cs.cliques_cons[clique_i])
            println(io, "           ", cs.cliques_cons[clique_i][con_j], " :  with $(length(cs.cliques_idcs_bases[clique_i][con_j])) basis monomials")
        end
    end
end

"""
    get_correlative_graph(ordered_vars::Vector{Variable}, obj::Polynomial{T}, cons::Vector{Polynomial{T}}, order::Int) where {T}

Constructs a correlative sparsity graph from polynomial optimization problem components.

# Arguments
- `ordered_vars::Vector{Variable}`: Variables in the order to appear in the graph
- `obj::Polynomial{T}`: Objective polynomial
- `cons::Vector{Polynomial{T}}`: Constraint polynomials
- `order::Int`: Order of the moment relaxation

# Returns
- `SimpleGraph`: Graph representing variable correlations
"""
function get_correlative_graph(ordered_vars::Vector{Variable}, obj::Polynomial{T}, cons::Vector{Polynomial{T}}, order::Int) where {T}
    # NOTE: code will be buggy is ordered_vars is not the same as the one reference in other functions
    @assert issorted(ordered_vars) "Variables must be sorted"

    nvars = length(ordered_vars)
    G = SimpleGraph(nvars)

    # find index of all unique variables in polynomial/monomial p
    vmap(p) = map(v -> findfirst(==(v), ordered_vars), unique!(variables(p)))

    map(mono -> add_clique!(G, vmap(mono)), obj.monos)

    for poly in cons
        # for clearer logic, I didn't combine the two branches
        if order == ceil(Int, maxdegree(poly) // 2)
            # if objective or order too large, each term forms a clique
            map(mono -> add_clique!(G, vmap(mono)), poly.monos)
        else
            # NOTE: if is constraint and order not too large, all variables in the constraint forms a clique
            # this ensures each "small" constraint is in a clique ?
            add_clique!(G, vmap(poly))
        end
    end
    return G
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

"""
    assign_constraint(cliques::Vector{Vector{Variable}}, cons::Vector{Polynomial{T}}) where {T}

Assigns constraints to cliques based on variable support.

# Arguments
- `cliques::Vector{Vector{Variable}}`: Variable cliques
- `cons::Vector{Polynomial{T}}`: Constraint polynomials

# Returns
- `Tuple{Vector{Vector{Int}}, Vector{Int}}`: Tuple containing:
  - Constraint indices for each clique
  - Global constraint indices not captured by any single clique
"""
function assign_constraint(cliques::Vector{Vector{Variable}}, cons::Vector{Polynomial{T}}) where {T}
    # assign each constraint to a clique
    # there might be constraints that are not captured by any single clique,
    # NOTE: we ignore this constraint. This should only occur at lower order of relaxation.

    # clique_cons: vector of vector of constraints index, each belong to a clique
    clique_cons = map(cliques) do clique
        findall(g -> issubset(unique!(variables(g)), clique), cons)
    end
    return clique_cons, setdiff(1:length(cons), union(clique_cons...))
end

"""
    correlative_sparsity(pop::PolyOpt{T}, order::Int, elim_algo::EliminationAlgorithm) where {T}

Decomposes a polynomial optimization problem into a correlative sparsity pattern by identifying
variable cliques and assigning constraints to cliques, enabling block-structured semidefinite relaxations.

# Arguments
- `pop::PolyOpt{T}`: Polynomial optimization problem containing objective, constraints, and variables
- `order::Int`: Order of the moment relaxation
- `elim_algo::EliminationAlgorithm`: Algorithm for clique tree elimination

# Returns
- `CorrelativeSparsity`: Structure containing:
  - `cliques`: Groups of variables that form cliques in the sparsity graph
  - `cliques_cons`: Constraint indices assigned to each clique
  - `global_cons`: Constraint indices not captured by any single clique
  - `cliques_idcs_bases`: Monomial bases for indexing moment/localizing matrices within each clique
"""
function correlative_sparsity(pop::PolyOpt{T}, order::Int, elim_algo::EliminationAlgorithm) where {T}
    cliques = map(x -> pop.variables[x], clique_decomp(get_correlative_graph(pop.variables, pop.objective, [pop.eq_constraints, pop.ineq_constraints], order), elim_algo))

    cliques_cons, global_cons = assign_constraint(cliques, pop.constraints)

    reduce_func = prod ∘ reducer(pop)
    # get the operators needed to index columns of moment/localizing mtx in each clique
    # depending on the clique's varaibles each is slightly different
    cliques_idx_basis = map(zip(cliques, cliques_cons)) do (clique, clique_cons)
        # get the basis of the moment matrix in a clique, then sort it
        [[sorted_unique(reduce_func.(get_basis(sort(clique, rev=true), order)))]; map(b -> sorted_unique(reduce_func.(b)), get_basis.(Ref(sort(clique, rev=true)), order .- ceil.(Int, maxdegree.(pop.constraints[clique_cons]) / 2)))]
    end

    return CorrelativeSparsity(cliques, cliques_cons, global_cons, cliques_idx_basis)
end


"""
    TermSparsity

Structure representing term sparsity information for polynomial optimization.

# Fields
- `term_sparse_graph_supp::Vector{Monomial}`: Support of the term sparsity graph
- `block_bases::Vector{Vector{Monomial}}`: Bases of moment/localizing matrices in each block
"""
struct TermSparsity
    term_sparse_graph_supp::Vector{Monomial}
    block_bases::Vector{Vector{Monomial}}
end

"""
    get_term_sparsity_graph(cons_support::Vector{Monomial}, activated_supp::Vector{Monomial}, basis::Vector{Monomial})

Constructs a term sparsity graph for polynomial constraints.

# Arguments
- `cons_support::Vector{Monomial}`: Support monomials of constraints
- `activated_supp::Vector{Monomial}`: Support from previous iterations
- `basis::Vector{Monomial}`: Basis used to index the moment matrix

# Returns
- `SimpleGraph`: Term sparsity graph
"""
function get_term_sparsity_graph(cons_support::Vector{Monomial}, activated_supp::Vector{Monomial}, basis::Vector{Monomial})
    nterms = length(basis)
    G = SimpleGraph(nterms)
    sorted_activated_supp = sort(activated_supp)
    for i in 1:nterms, j in i+1:nterms
        for supp in cons_support
            if symmetric_canonicalize(neat_dot(basis[i], supp * basis[j])) in sorted_activated_supp
                add_edge!(G, i, j)
                continue
            end
        end
    end
    return G
end

"""
    iterate_term_sparse_supp(activated_supp::Vector{Monomial}, poly::Polynomial, basis::Vector{Monomial}, elim_algo::EliminationAlgorithm)

Iteratively computes term sparsity support for a polynomial.

# Arguments
- `activated_supp::Vector{Monomial}`: Currently activated support monomials
- `poly::Polynomial`: Input polynomial
- `basis::Vector{Monomial}`: Basis monomials
- `elim_algo::EliminationAlgorithm`: Elimination algorithm for clique decomposition

# Returns
- `TermSparsity`: Term sparsity structure containing graph support and block bases
"""
function iterate_term_sparse_supp(activated_supp::Vector{Monomial}, poly::Polynomial, basis::Vector{Monomial}, elim_algo::EliminationAlgorithm)
    F = get_term_sparsity_graph(poly.monos, activated_supp, basis)
    if !(elim_algo isa AsIsElimination)
        blocks = clique_decomp(F, elim_algo)
        map(block -> add_clique!(F, block), blocks)
    else
        blocks = connected_components(F)
    end
    return TermSparsity(term_sparsity_graph_supp(F, basis, poly), map(x -> basis[x], blocks))
end

"""
    term_sparsity_graph_supp(G::SimpleGraph, basis::Vector{Monomial}, g::Polynomial)

Computes the support of a term sparsity graph for a given polynomial.

# Arguments
- `G::SimpleGraph`: Term sparsity graph
- `basis::Vector{Monomial}`: Basis monomials
- `g::Polynomial`: Input polynomial

# Returns
- `Vector{Monomial}`: Support monomials for the term sparsity graph
"""
function term_sparsity_graph_supp(G::SimpleGraph, basis::Vector{Monomial}, g::Polynomial)
    # following (10.4) in Sparse Polynomial Optimization: Theory and Practise
    # NOTE: Do I need to symmetric canonicalize it?
    # TODO: add reduce! here
    gsupp(a, b) = map(g_supp -> neat_dot(a, g_supp * b), g.monos)
    return union([gsupp(basis[v], basis[v]) for v in vertices(G)]..., [gsupp(basis[e.src], basis[e.dst]) for e in edges(G)]...)
end
