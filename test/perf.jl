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
    cs_algo = MF()
    ts_algo = MMD()
    solver_config = SolverConfig(optimizer=Mosek.Optimizer; mom_order=3, cs_algo=cs_algo, ts_algo=ts_algo)

    mom_order = iszero(solver_config.mom_order) ? maximum([ceil(Int, maxdegree(poly) / 2) for poly in [pop.objective; pop.constraints]]) : solver_config.mom_order

    using NCTSSoS: correlative_sparsity
    @benchmark corr_sparsity = correlative_sparsity($pop, $mom_order, $cs_algo)
    corr_sparsity = correlative_sparsity(pop, mom_order, solver_config.cs_algo)

    using DynamicPolynomials: coefficients, monomials, effective_variables

    #  823.374 μs (22103 allocations: 2.07 MiB)
    cliques_objective = [reduce(+, [issubset(effective_variables(mono), clique) ? coef * mono : zero(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

    using NCTSSoS: sorted_union
    using DynamicPolynomials: Monomial
    # 1.359 ms (60361 allocations: 3.70 MiB)
    initial_activated_supp = [sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, pop.constraints[cons_idx]; init=eltype(monomials(obj_part))[]), [neat_dot(b, b) for b in idcs_bases[1]])
                              for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)]

    using NCTSSoS: iterate_term_sparse_supp
    # 4.211 s (118776605 allocations: 8.11 GiB)
    cliques_term_sparsities = map(zip(initial_activated_supp, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
        [iterate_term_sparse_supp(activated_supp, poly, basis, solver_config.ts_algo) for (poly, basis) in zip([one(pop.objective); pop.constraints[cons_idx]], idcs_bases)]
    end

    using NCTSSoS: moment_relax
    # 69.350 ms (3497593 allocations: 154.94 MiB)
    moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, corr_sparsity.global_cons, cliques_term_sparsities)

    using NCTSSoS: sos_dualize
    #   124.614 ms (5812129 allocations: 264.22 MiB)
    @btime begin
    sos_problem = sos_dualize(moment_problem)
    end
    set_optimizer(sos_problem.model, solver_config.optimizer)
end