struct PolyOptResult{T,P,M}
    objective::T # support for high precision solution
    corr_sparsity::CorrelativeSparsity{P,M}
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
    # should contain objective and moment matrix or gram matrix for manually checking what happened
end

function Base.show(io::IO, result::PolyOptResult)
    println(io, "Objective: ", result.objective)
    show(io, result.corr_sparsity)
    println(io, "Term Sparsity:")
    for (i, sparsities) in enumerate(result.cliques_term_sparsities)
        println(io, "Clique $i:")
        println(io, "   Moment Matrix:")
        println(io, sparsities[1])
        println(io, "   Localizing Matrix:")
        for sparsity in sparsities[2:end]
            show(io, sparsity)
        end
    end
end

"""
    SolverConfig(; optimizer, mom_order, cs_algo=NoElimination(), ts_algo=NoElimination())

Configuration for solving polynomial optimization problems.

# Keyword Arguments
- `optimizer` (required): The optimizer to use for solving the SDP problem (e.g. Clarabel.Optimizer)
- `mom_order::Int`: The order of the moment relaxation (default: 0)
- `cs_algo::EliminationAlgorithm`: Algorithm for correlative sparsity exploitation (default: NoElimination())
- `ts_algo::EliminationAlgorithm`: Algorithm for term sparsity exploitation (default: NoElimination())


# Examples
```jldoctest; setup=:(using NCTSSoS, Clarabel)
julia> solver_config = SolverConfig(optimizer=Clarabel.Optimizer, mom_order=2) # default elimination algorithms
SolverConfig(Clarabel.MOIwrapper.Optimizer, 2, NoElimination(), NoElimination())
```

"""
@kwdef struct SolverConfig
    optimizer
    mom_order::Int = 0
    cs_algo::EliminationAlgorithm = NoElimination()
    ts_algo::EliminationAlgorithm = NoElimination()
end

# TODO: add user interface
# learn from OMEinsum.jl
# consider adding Solver Interface
# consider obtaining enough information on Moment matrix etc to check if problem solved correctly
# prev_ans::Union{Nothing,PolyOptResult{C,T}}=nothing
function cs_nctssos(pop::PolyOpt{T}, solver_config::SolverConfig; dualize::Bool=true) where {T}
    sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, is_projective=pop.is_projective, is_unipotent=pop.is_unipotent)
    mom_order = iszero(solver_config.mom_order) ? maximum([ceil(Int, maxdegree(poly) / 2) for poly in [pop.objective; pop.eq_constraints; pop.ineq_constraints]]) : solver_config.mom_order

    corr_sparsity = correlative_sparsity(pop, mom_order, solver_config.cs_algo)

    cliques_objective = [reduce(+, [issubset(sort!(variables(mono)), clique) ? coef * mono : zero(mono) for (coef, mono) in zip(pop.objective.coeffs, pop.objective.monos)]) for clique in corr_sparsity.cliques]

    initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base, sa)
    end

    cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo, sa)
    end

    moment_problem = moment_relax(pop, corr_sparsity, cliques_term_sparsities)
    if dualize
        sos_problem = sos_dualize(moment_problem)
        set_optimizer(sos_problem.model, solver_config.optimizer)
        optimize!(sos_problem.model)
        return PolyOptResult(objective_value(sos_problem.model), corr_sparsity, cliques_term_sparsities)
    else
        set_optimizer(moment_problem.model, solver_config.optimizer)
        optimize!(moment_problem.model)
        return PolyOptResult(objective_value(moment_problem.model), corr_sparsity, cliques_term_sparsities)
    end
end

function cs_nctssos_higher(pop::PolyOpt{T}, prev_res::PolyOptResult, solver_config::SolverConfig; dualize::Bool=true) where {T}
    sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, is_unipotent=pop.is_unipotent, is_projective=pop.is_projective)
    initial_activated_supps = [sorted_union([poly_term_sparsity.term_sparse_graph_supp for poly_term_sparsity in term_sparsities]...)
                               for term_sparsities in prev_res.cliques_term_sparsities]

    prev_corr_sparsity = prev_res.corr_sparsity
    cliques_term_sparsities = map(zip(initial_activated_supps, prev_corr_sparsity.clq_cons, prev_corr_sparsity.clq_mom_mtx_bases, prev_corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, prev_corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo, sa)
    end


    moment_problem = moment_relax(pop, prev_res.corr_sparsity, cliques_term_sparsities)
    if !dualize
        set_optimizer(moment_problem.model, solver_config.optimizer)
        optimize!(moment_problem.model)
        return PolyOptResult(objective_value(moment_problem.model), prev_res.corr_sparsity, cliques_term_sparsities)
    else
        sos_problem = sos_dualize(moment_problem)
        set_optimizer(sos_problem.model, solver_config.optimizer)
        optimize!(sos_problem.model)
        return PolyOptResult(objective_value(sos_problem.model), prev_res.corr_sparsity, cliques_term_sparsities)
    end
end


# TODO: add user interface
# learn from OMEinsum.jl
# consider adding Solver Interface
# consider obtaining enough information on Moment matrix etc to check if problem solved correctly
# prev_ans::Union{Nothing,PolyOptResult{C,T}}=nothing
function cs_nctssos(spop::PolyOpt{NCStatePolynomial{T},OBJ}, solver_config::SolverConfig; dualize::Bool=true) where {T,OBJ}

    mom_order = iszero(solver_config.mom_order) ? maximum([ceil(Int, maxdegree(poly) / 2) for poly in [spop.objective; spop.eq_constraints; spop.ineq_constraints]]) : solver_config.mom_order

    reducer_func = reducer(spop)

    corr_sparsity = correlative_sparsity(spop, mom_order, solver_config.cs_algo)

    cliques_objective = [reduce(+, [issubset(sort!(variables(t[2])), clique) ? t[1] * t[2] : zero(t[2]) for t in terms(spop.objective)]) for clique in corr_sparsity.cliques]

    initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base, reducer_func)
    end

    cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo, reducer_func)
    end

    moment_problem = moment_relax(spop, corr_sparsity, cliques_term_sparsities)
    if dualize
        sos_problem = sos_dualize(moment_problem)
        set_optimizer(sos_problem.model, solver_config.optimizer)
        optimize!(sos_problem.model)
        return PolyOptResult(objective_value(sos_problem.model), corr_sparsity, cliques_term_sparsities)
    else
        set_optimizer(moment_problem.model, solver_config.optimizer)
        optimize!(moment_problem.model)
        return PolyOptResult(objective_value(moment_problem.model), corr_sparsity, cliques_term_sparsities)
    end
end
