using NCTSSoS, NCTSSoS.FastPolynomials
using MosekTools
using JuMP

function cs_nctssos_with_entry(pop::OP, solver_config::SolverConfig, entry_constraints::Vector{Polynomial{T}}; dualize::Bool=true) where {T,P<:Polynomial{T},OP<:NCTSSoS.OptimizationProblem{P}}

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

N = 10 
T = ComplexF64

# this is supposed to be upper bounds because it's from DMRG
energy_lower_bounds = [
     -0.43286676510889616,
     -0.39861081639315044,
     -0.3750000000000001,
     -0.403134244708691,
     -0.46685450336963996,
     -0.5358985722626579,
     -0.6067542029652743,
     -0.678582898886551,
     -0.7510320744136272,
     -0.8239121834672108,
     -0.8971075269023185
]

energy_upper_bounds = [


]

for i in 1:11
     J1, J2 = 1.0, 0.1 + (i - 1) * 0.2

     @ncpolyvar x[1:N] y[1:N] z[1:N]

     obj = one(T) * x[1] * x[2]

     ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] + T(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

     eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

     solver_config = SolverConfig(optimizer=SOLVER, order=2)

     pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

     res = cs_nctssos(pop, solver_config)

     @assert is_solved_and_feasible(res.model)
     push!(energy_upper_bounds, res.objective)
end

energy_upper_bounds

pop_l = cpolyopt(obj; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

pop_u = cpolyopt(-obj; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=SOLVER, order=3)

# would it be better to use (H - (up+lb)/2)^2 = sqrt(up-lb) ?
ineq_cons = [ham - energy_lower_bounds[i] * N, energy_upper_bounds[i] * N - ham]
# ineq_cons = [-(ham - (energy_lower_bounds[i] * N + energy_upper_bounds[i] * N) / 2)^2 + sqrt(energy_upper_bounds[i] - energy_lower_bounds[i])]

res_l = cs_nctssos_with_entry(pop_l, solver_config, ineq_cons; dualize=true)

res_u = cs_nctssos_with_entry(pop_u, solver_config, ineq_cons; dualize=true)

res_l.objective/4
-res_u.objective/4

termination_status(res_l.model)
is_solved_and_feasible(res_l.model)

termination_status(res_u.model)
is_solved_and_feasible(res_u.model)