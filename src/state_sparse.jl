# cliques: grouping of variables into sets of variables
# cliques_cons: groups constraints according to cliques,
# constraints in a clique only has support on corresponding variables
# discarded_cons: constraints that are not in any clique
# cliques_idcs_bases: within each clique, the vectors of monomials used to index moment/localizing matrices
# TODO: merge this with CorrelativeSparsity by adding a parameter
struct StateCorrelativeSparsity{V,M}
    cliques::Vector{Vector{Variable{V,M}}}
    cliques_cons::Vector{Vector{Int}}
    # FIXME: add test case for difference
    global_cons::Vector{Int}
    cliques_idcs_bases::Vector{Vector{Vector{NCStateWord{V,M}}}}
end

# ordered_vars: variables in the order to be appeared in graph
# polys: objective + constraints, order is important
# order: order of the moment problem
# TODO: this can also be merged with get_correlative_graph
function get_correlative_graph(ordered_vars::Vector{Variable{V,M}}, obj::StatePolynomialOp{V,M,T}, cons::Vector{StatePolynomialOp{V,M,T}}, order::Int) where {V,M,T}
    # NOTE: code will be buggy is ordered_vars is not the same as the one reference in other functions
    # @assert issorted(ordered_vars, rev=true) "Variables must be sorted"

    nvars = length(ordered_vars)
    G = SimpleGraph(nvars)

    # find index of all unique variables in polynomial/monomial p
    vmap(p) = map(v -> findfirst(==(v), ordered_vars), unique!(effective_variables(p)))

    map(mono -> add_clique!(G, vmap(mono)), monomials(obj))

    for poly in cons
        # for clearer logic, I didn't combine the two branches
        if order == ceil(Int, maxdegree(poly) // 2)
            # if objective or order too large, each term forms a clique
            map(mono -> add_clique!(G, vmap(mono)), monomials(poly))
        else
            # NOTE: if is constraint and order not too large, all variables in the constraint forms a clique
            # this ensures each "small" constraint is in a clique ?
            add_clique!(G, vmap(poly))
        end
    end
    return G
end

function assign_constraint(cliques::Vector{Vector{Variable{V,M}}}, cons::Vector{StatePolynomialOp{V,M,T}}) where {V,M,T}
    # assign each constraint to a clique
    # there might be constraints that are not captured by any single clique,
    # NOTE: we ignore this constraint. This should only occur at lower order of relaxation.

    # clique_cons: vector of vector of constraints index, each belong to a clique
    clique_cons = map(cliques) do clique
        findall(g -> issubset(unique!(effective_variables(g)), clique), cons)
    end
    return clique_cons, setdiff(1:length(cons), union(clique_cons...))
end

function correlative_sparsity(pop::StatePolyOpt{V,M,T}, order::Int, elim_algo::EliminationAlgorithm) where {V,M,T}
    cliques = map(x -> pop.variables[x], clique_decomp(get_correlative_graph(pop.variables, pop.objective, pop.constraints, order), elim_algo))

    cliques_cons, global_cons = assign_constraint(cliques, pop.constraints)

    reduce_func = reducer(pop)
    # get the operators needed to index columns of moment/localizing mtx in each clique
    # depending on the clique's varaibles each is slightly different
    cliques_idx_basis = map(zip(cliques, cliques_cons)) do (clique, clique_cons)
        # get the basis of the moment matrix in a clique, then sort it
        [[sorted_unique(reduce_func.(get_state_basis(sort(clique, rev=true), order, reduce_func)))]; map(b -> sorted_unique(reduce_func.(b)), get_state_basis.(Ref(sort(clique, rev=true)), order .- ceil.(Int, maxdegree.(pop.constraints[clique_cons]) / 2)), Ref(reduce_func))]
    end

    return StateCorrelativeSparsity{V,M}(cliques, cliques_cons, global_cons, cliques_idx_basis)
end

# term_sparse_graph_supp: support of the current term sparsity graph for an obj/cons
# block_bases: the bases of the moment/localizing matrix in each clique of term sparse graph
struct StateTermSparsity{V,M}
    term_sparse_graph_supp::Vector{NCStateWord{V,M}}
    block_bases::Vector{Vector{NCStateWord{V,M}}}
end

# porting nccpop.jl's  get_graph
# constructs the graph according to (7.5) and (7.14) together
# activated_supp: support of objective, constraint and their corresponding term sparsity graph in previous iteration (7.14)
# basis: basis used to index the moment matrix
function get_term_sparsity_graph(cons_support::Vector{NCStateWord{V,M}}, activated_supp::Vector{NCStateWord{V,M}}, basis::Vector{NCStateWord{V,M}}) where {V,M}
    nterms = length(basis)
    G = SimpleGraph(nterms)
    as = expval.(activated_supp)
    for i in 1:nterms, j in i+1:nterms
        for supp in cons_support
            # interm = symmetric_canonicalize(neat_dot(basis[i], supp * basis[j]))
            interm = neat_dot(basis[i], supp*basis[j])
            if cyclic_canonicalize(expval(interm)) in as
                add_edge!(G, i, j)
                continue
            end
        end
    end
    return G
end

# returns: F (the chordal graph), blocks in basis
function iterate_term_sparse_supp(activated_supp::Vector{NCStateWord{V,M}}, poly::StatePolynomialOp, basis::Vector{NCStateWord{V,M}}, elim_algo::EliminationAlgorithm) where {V,M}
    F = get_term_sparsity_graph(collect(monomials(poly)), activated_supp, basis)
    if !(elim_algo isa AsIsElimination)
        blocks = clique_decomp(F, elim_algo)
        map(block -> add_clique!(F, block), blocks)
    else
        blocks = connected_components(F)
    end
    return StateTermSparsity(term_sparsity_graph_supp(F, basis, poly), map(x -> basis[x], blocks))
end

# supp(G,g): monomials that are either v^† g_supp v where v is a vertex in G, or β^† g_supp γ where {β,γ} is an edge in G following (10,4)
# given term sparsity graph G, which terms needs to be considered as a variable for describing the localizing/moment matrix with respect to g
function term_sparsity_graph_supp(G::SimpleGraph, basis::Vector{NCStateWord{V,M}}, g::StatePolynomialOp) where {V,M}
    # following (10.4) in Sparse Polynomial Optimization: Theory and Practise
    # NOTE: Do I need to symmetric canonicalize it?
    # TODO: add reduce! here
    gsupp(a, b) = map(g_supp -> neat_dot(a, g_supp * b), monomials(g))
    return union([gsupp(basis[v], basis[v]) for v in vertices(G)]..., [gsupp(basis[e.src], basis[e.dst]) for e in edges(G)]...)
end
