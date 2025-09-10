using NCTSSoS, NCTSSoS.FastPolynomials
using MosekTools
using JuMP

 -0.11892547876100075
 -0.1999999999999994
 -0.288675134594812
 -0.38134355029701733
 -0.47627352424814867
 -0.5725815626252594
 -0.6697825646477007
 -0.7675918792439959
 -0.8658328118479361
 -0.9643904066487076

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

N = 3
T = ComplexF64
i = 1 

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

# energy_lower_bounds = [
# -0.44667162694019064,
# -0.40833333333333327,
# -0.37499999999999994,
# -0.42500000000000004,
# -0.47499999999999987,
# -0.5249999999999999,
# -0.5749999999999998,
# -0.6249999999999997,
# -0.6749999999999998,
# -0.7249999999999998,
# -0.7749999999999998
# ]

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

# energy_upper_bounds = [
#      -0.4466716248248655,
#      -0.40833333318460446,
#      -0.3749999997829185,
#      -0.4249999999582507,
#      -0.47499999962478495,
#      -0.5249999997928769,
#      -0.5749999999996329,
#      -0.6249999999996749,
#      -0.6749999999905304,
#      -0.7249999971735893,
#      -0.774999997626372
# ]

i = 1
J1, J2 = 1.0, 0.1 + (i - 1) * 0.2

@ncpolyvar x[1:N] y[1:N] z[1:N]

obj = one(T) * x[1] * x[2]

ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [x, y, z] for i in 1:N) + sum(T(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in [1])

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])


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

using Yao
using LinearAlgebra

N = 3
J1, J2 = 1.0, 0.1

ham = sum(J1 / 4 * kron(N, i => op, mod1(i + 1, N) => op) for op in (X, Y, Z) for i in 1:N) + sum(J2 / 4 * kron(N, i => op, mod1(i + 2, N) => op) for op in (X, Y, Z) for i in [1])

s1s2 = Matrix(kron(N,1=>X,2=>X))/4

evals, eigvecs = eigen(Matrix(ham))

# has degeneracy
evals

evals ./ N
evals[end]

s1s2s = map(1:length(evals)) do which_state
    # eigvecs[:, which_state]' * Matrix(ham) * eigvecs[:, which_state] / N
    real(eigvecs[:, which_state]' * s1s2 * eigvecs[:, which_state])
end
