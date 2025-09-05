using NCTSSoS, NCTSSoS.FastPolynomials
using Clarabel
# using MosekTools
using JuMP

SOLVER = optimizer_with_attributes(Clarabel.Optimizer, "direct_solve_method" => :cudss)

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

# SOLVER = optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-8, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => 1e-8, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1e-8, "MSK_IPAR_NUM_THREADS" => 0)

N = 4
T = ComplexF64

# this is supposed to be upper bounds because it's from DMRG
energy_lower_bounds = [
     -0.475,
     -0.42500000000000004,
     -0.3749999999999999,
     -0.5249999999999996,
     -0.6749999999999997,
     -0.8250000000000001,
     -0.9749999999999995,
     -1.125,
     -1.2749999999999995,
     -1.4249999999999994,
     -1.5749999999999997]

energy_upper_bounds = [
     -0.47499999725973996,
     -0.424999998456548,
     -0.3749999999529608,
     -0.5249999999094307,
     -0.674999999842614,
     -0.8249999985406258,
     -0.9749999998616854,
     -1.1249999996288769,
     -1.2749999993449743,
     -1.4249999999641514,
     -1.574999999996572
]

lower_bounds = Float64[]
upper_bounds = Float64[]

# for i in 1:11
i = 1
J1, J2 = 1.0, 0.1 + (i - 1) * 0.2

@ncpolyvar x[1:N] y[1:N] z[1:N]

obj = one(T) * x[1] * x[2]

ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] + T(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

ineq_cons = [ham - energy_lower_bounds[i] * N]

pop_l = cpolyopt(obj; eq_constraints=eq_cons, ineq_constraints=ineq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

pop_u = cpolyopt(-obj; eq_constraints=eq_cons, ineq_constraints=ineq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=SOLVER, order=2)

res = cs_nctssos(pop_l, solver_config)

# ineq_cons = [energy_upper_bounds[i] * N - ham]

# res_l = cs_nctssos_with_entry(pop_l, solver_config, ineq_cons; dualize=true)
# res_u = cs_nctssos_with_entry(pop_u, solver_config, ineq_cons; dualize=true)

is_solved_and_feasible(res_l.model) || error("Lower bound problem not solved or feasible at J2 = $J2")
is_solved_and_feasible(res_u.model) || error("Upper bound problem not solved or feasible at J2 = $J2")

push!(lower_bounds, res_l.objective / 4)
push!(upper_bounds, -res_u.objective / 4)
# end

# the larger the J2 the harder it is to obtain good bounds of s0s1
# at J2 = 0.9, this formulation (order = 2) fails to bound s0s1

using DelimitedFiles

open("res.txt","w") do io
     writedlm(io, hcat(lower_bounds, upper_bounds))
end
