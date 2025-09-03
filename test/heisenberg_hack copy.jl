using NCTSSoS, NCTSSoS.FastPolynomials
using MosekTools
using JuMP

# hacky way of adding entry-wise inequality
function cs_nctssos_with_entry!(pop::OP, solver_config::SolverConfig, entry_constraints::Vector{Polynomial{T}}; dualize::Bool=true) where {T,P<:Polynomial{T},OP<:NCTSSoS.OptimizationProblem{P}}

   sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, is_projective=pop.is_projective, is_unipotent=pop.is_unipotent)
   order = iszero(solver_config.order) ? maximum([ceil(Int, maxdegree(poly) / 2) for poly in [pop.objective; pop.eq_constraints; pop.ineq_constraints]]) : solver_config.order

   corr_sparsity = NCTSSoS.correlative_sparsity(pop, order, solver_config.cs_algo)

     cliques_objective = [reduce(+, [issubset(sort!(variables(mono)), clique) ? coef * mono : zero(coef) * one(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

   initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        NCTSSoS.init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base, sa)
   end

   cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        NCTSSoS.term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo, sa)
   end

   moment_problem = NCTSSoS.moment_relax(pop, corr_sparsity, cliques_term_sparsities)

   for c in entry_constraints
     push!(moment_problem.constraints,(:HPSD, [c;;]))
   end

   (pop isa NCTSSoS.ComplexPolyOpt{P} && !dualize) && error("Solving Moment Problem for Complex Poly Opt is not supported")
   problem_to_solve = !dualize ? moment_problem : NCTSSoS.sos_dualize(moment_problem)

   set_optimizer(problem_to_solve.model, solver_config.optimizer)
   optimize!(problem_to_solve.model)
   return NCTSSoS.PolyOptResult(objective_value(problem_to_solve.model), corr_sparsity, cliques_term_sparsities, problem_to_solve.model)
end

SOLVER = optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-8, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => 1e-8, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1e-8, "MSK_IPAR_NUM_THREADS" => 0)

N = 1
T = ComplexF64

@ncpolyvar x[1:N] y[1:N] z[1:N]

solver_config = SolverConfig(optimizer=SOLVER, order=4)

op = one(T) * x[1]

ham = one(T) * z[1]

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

lower_bound_energy = -1.000000000001
upper_bound_energy = -0.999999999999
ineq_cons = [ham - lower_bound_energy, upper_bound_energy - ham]

pop = cpolyopt(op; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)
res = cs_nctssos_with_entry!(pop, solver_config, ineq_cons; dualize=true)

@assert isapprox(res.objective, 0; atol=1e-5)
termination_status(res.model)

@assert is_solved_and_feasible(res.model)

for lower_bound_energy in 1.0:-0.1:-1.0
     ineq_cons = [ham - lower_bound_energy]

     pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

     res = cs_nctssos_with_entry!(pop, solver_config, ineq_cons; dualize=true)

     @assert isapprox(res.objective ,lower_bound_energy;atol=1e-5)

     @show termination_status(res.model)

     @assert is_solved_and_feasible(res.model)
end