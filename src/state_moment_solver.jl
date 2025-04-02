# V: is the varaibles commuting
# M: Ordering of variables
# T: type of the coefficients
# monomap: map from monomials in DynamicPolynomials to variables in JuMP
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
function moment_relax(pop::StatePolyOpt{V,M,T}, mom_order::Int) where {V,M,T}
    # NOTE: objective and constraints may have integer coefficients, but popular JuMP solvers does not support integer coefficients
    # left type here to support BigFloat model for higher precision
    model = GenericModel{T}()

    reduce_func = reducer(pop)

    # the union of clique_total_basis
    column_basis = get_state_basis(pop.variables, mom_order)
    total_basis = sorted_unique(vec(
        [
        expval(reduce_func(neat_dot(a, b)))
        for a in column_basis, b in column_basis]
    ))

    # # map the monomials to JuMP variables, the first variable must be 1
    @variable(model, y[1:length(total_basis)])
    @constraint(model, first(y) == 1)
    monomap = Dict(zip(total_basis, y))

    constraint_matrices = [mapreduce(vcat, zip([false; pop.is_equality], [one(pop.objective); pop.constraints])) do (iseq, cons)
        constrain_moment_matrix!(model, cons, get_state_basis(pop.variables, mom_order - fld(degree(cons), 2)), monomap, iseq ? Zeros() : PSDCone(), reduce_func)
    end]

    @objective(model, Min, substitute_variables(mapreduce(p -> p.coef * reduce_func(p.ncstate_word), +, terms(pop.objective); init=zero(pop.objective)), monomap))

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
        substitute_variables(sum([poly_term.coef * reduce_func(neat_dot(row_idx, poly_term.ncstate_word * col_idx)) for poly_term in terms(poly)];init=zero(poly)), monomap) for row_idx in local_basis, col_idx in local_basis
    ]
    return @constraint(model, moment_mtx in cone)
end
