using BenchmarkTools
using MosekTools, NCTSSoS

@testset "Documentation Example" begin
    n = 5
    @ncpolyvar x[1:n]
    f = 0.0
    for i = 1:n
        jset = max(1, i - 5):min(n, i + 1)
        jset = setdiff(jset, i)
        f += (2x[i] + 5 * x[i]^3 + 1)^2
        f -= sum([4x[i] * x[j] + 10x[i]^3 * x[j] + 2x[j] + 4x[i] * x[j]^2 + 10x[i]^3 * x[j]^2 + 2x[j]^2 for j in jset])
        f += sum([x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset])
    end

	cons = typeof(f)[]
    for i = 1:n
        push!(cons, 1 - x[i]^2)
        push!(cons, x[i] - 1 / 3)
    end

    pop =  PolyOpt(f; constraints = cons);
    solver_config = SolverConfig(optimizer=Mosek.Optimizer; mom_order=3, cs_algo=MF(), ts_algo=MMD())

    @benchmark result_cs_ts = cs_nctssos($pop, $solver_config)
    result_cs_ts = cs_nctssos(pop, solver_config)

    @benchmark result_cs_ts_higher = cs_nctssos_higher($pop, $result_cs_ts, $solver_config)
end

@testset "Components" begin
    n = 5
    @ncpolyvar x[1:n]
    f = 0.0
    for i = 1:n
        jset = max(1, i - 5):min(n, i + 1)
        jset = setdiff(jset, i)
        f += (2x[i] + 5 * x[i]^3 + 1)^2
        f -= sum([4x[i] * x[j] + 10x[i]^3 * x[j] + 2x[j] + 4x[i] * x[j]^2 + 10x[i]^3 * x[j]^2 + 2x[j]^2 for j in jset])
        f += sum([x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset])
    end

	cons = typeof(f)[]
    for i = 1:n
        push!(cons, 1 - x[i]^2)
        push!(cons, x[i] - 1 / 3)
    end

    pop =  PolyOpt(f; constraints = cons);
    solver_config = SolverConfig(optimizer=Mosek.Optimizer; mom_order=3, cs_algo=MF(), ts_algo=MMD())

    mom_order = iszero(solver_config.mom_order) ? maximum([ceil(Int, maxdegree(poly) / 2) for poly in [pop.objective; pop.constraints]]) : solver_config.mom_order

    using NCTSSoS: correlative_sparsity
    @btime begin
    corr_sparsity = correlative_sparsity(pop, mom_order, solver_config.cs_algo)
    end

    cliques_objective = [reduce(+, [issubset(effective_variables(mono), clique) ? coef * mono : zero(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

    # prepare the support for each term sparse localizing moment
    initial_activated_supp = [sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, pop.constraints[cons_idx]; init=Monomial{V,M}[]), [neat_dot(b, b) for b in idcs_bases[1]])
                              for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)]

    cliques_term_sparsities = map(zip(initial_activated_supp, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
        [iterate_term_sparse_supp(activated_supp, poly, basis, solver_config.ts_algo) for (poly, basis) in zip([one(pop.objective); pop.constraints[cons_idx]], idcs_bases)]
    end

    moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, corr_sparsity.global_cons, cliques_term_sparsities)
    sos_problem = sos_dualize(moment_problem)
    set_optimizer(sos_problem.model, solver_config.optimizer)


end