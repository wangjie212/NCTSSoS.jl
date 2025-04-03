# V: is the varaibles commuting
# M: Ordering of variables
# T: type of the coefficients
# monomap: map from monomials in DynamicPolynomials to variables in JuMP
# TODO: move StateWord{V,M} to parameter of MomentProblem
struct StateMomentProblem{V,M,T,CR<:ConstraintRef} <: OptimizationProblem
    model::GenericModel{T}
    constraints::Vector{CR}
    monomap::Dict{StateWord{V,M},GenericVariableRef{T}}  # TODO: maybe refactor.
    reduce_func::Function
end
function substitute_variables(poly::StatePolynomialOp{V,M,T}, wordmap::Dict{StateWord{V,M},GenericVariableRef{T}}) where {V,M,T}
    mapreduce(x -> (x.coef * wordmap[expval(x.ncstate_word)]), +, terms(poly))
end

# cliques_cons: groups constraints according to cliques,
# global_cons: constraints that are not in any single clique
# cliques_term_sparsities: each clique, each obj/constraint, each ts_clique, each basis needed to index moment matrix
# FIXME: should I use CorrelativeSparsity here instead of cliques_cons and global_cons
function moment_relax(pop::StatePolyOpt{V,M,T}, cliques_cons::Vector{Vector{Int}}, global_cons::Vector{Int}, cliques_term_sparsities::Vector{Vector{StateTermSparsity{V,M}}}) where {V,M,T}
    # NOTE: objective and constraints may have integer coefficients, but popular JuMP solvers does not support integer coefficients
    # left type here to support BigFloat model for higher precision
    model = GenericModel{T}()

    reduce_func = reducer(pop)

    # the union of clique_total_basis
    total_basis = sorted_union(map(zip(cliques_cons, cliques_term_sparsities)) do (cons_idx, term_sparsities)
        union(vec(reduce(vcat, [
            map(monomials(poly)) do m
                expval(reduce_func(neat_dot(rol_idx, m * col_idx)))
            end
            for (poly, term_sparsity) in zip([one(pop.objective); pop.constraints[cons_idx]], term_sparsities) for basis in term_sparsity.block_bases for rol_idx in basis for col_idx in basis
        ])))
    end...)

    # # map the monomials to JuMP variables, the first variable must be 1
    @variable(model, y[1:length(total_basis)])
    @constraint(model, first(y) == 1)
    monomap = Dict(zip(total_basis, y))


    constraint_matrices =
        [mapreduce(vcat, zip(cliques_term_sparsities, cliques_cons)) do (term_sparsities, cons_idx)
                mapreduce(vcat, zip(term_sparsities, [one(pop.objective), pop.constraints[cons_idx]...], [false; pop.is_equality[cons_idx]])) do (term_sparsity, poly, is_eq)
                    map(term_sparsity.block_bases) do ts_sub_basis
                        constrain_moment_matrix!(
                            model,
                            poly,
                            ts_sub_basis,
                            monomap,
                            is_eq ? Zeros() : PSDCone(), reduce_func)
                    end
                end
            end
            map(global_cons) do global_con
                constrain_moment_matrix!(
                    model,
                    pop.constraints[global_con],
                    [one(pop.objective)],
                    monomap,
                    pop.is_equality[global_con] ? Zeros() : PSDCone(),
                    reduce_func
                )
            end]

    # constraint_matrices = [mapreduce(vcat, zip([false; pop.is_equality], [one(pop.objective); pop.constraints])) do (iseq, cons)
    #     constrain_moment_matrix!(model, cons, get_state_basis(pop.variables, mom_order - fld(degree(cons), 2)), monomap, iseq ? Zeros() : PSDCone(), reduce_func)
    # end]

    @objective(model, Min, substitute_variables(mapreduce(p -> p.coef * reduce_func(p.ncstate_word), +, terms(symmetric_canonicalize(pop.objective)); init=zero(pop.objective)), monomap))


    return StateMomentProblem(model, constraint_matrices, monomap, reduce_func)
end

function constrain_moment_matrix!(
    model::GenericModel{T},
    poly::StatePolynomialOp{V,M,T},
    local_basis::Vector{NCStateWord{V,M}},
    monomap::Dict{StateWord{V,M},GenericVariableRef{T}},
    cone, # FIXME: which type should I use?
    reduce_func::Function
) where {V,M,T}
    moment_mtx = [
        substitute_variables(sum([poly_term.coef * reduce_func(neat_dot(row_idx, poly_term.ncstate_word * col_idx)) for poly_term in terms(poly)]; init=zero(poly)), monomap) for row_idx in local_basis, col_idx in local_basis
    ]
    return @constraint(model, moment_mtx in cone)
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
    for i in 1:nterms, j in i+1:nterms
        for supp in cons_support
            if symmetric_canonicalize(neat_dot(basis[i], supp * basis[j])) in activated_supp
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
    blocks = clique_decomp(F, elim_algo)
    map(block -> add_clique!(F, block), blocks)
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
