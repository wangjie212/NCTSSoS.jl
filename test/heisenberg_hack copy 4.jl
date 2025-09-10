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

N = 3
T = ComplexF64

# this is supposed to be upper bounds because it's from DMRG
energy_lower_bounds = [
     -0.11892547876100075,
     -0.1999999999999994,
     -0.288675134594812,
     -0.38134355029701733,
     -0.47627352424814867,
     -0.5725815626252594,
     -0.6697825646477007,
     -0.7675918792439959,
     -0.8658328118479361,
     -0.9643904066487076]

energy_upper_bounds = Float64[
     -0.11892547875898078,
     -0.19999999979682634,
     -0.2886751335933035,
     -0.3813435496845692,
     -0.47627352422295716,
     -0.5725815625410057,
     -0.6697825646314844,
     -0.7675918790994536,
     -0.8658328114858912,
     -0.9643904052079089
]

s1s2vals = [
     -0.08130486019121048,
     -0.07142857142857142,
     -0.06100423396407306,
     -0.05235872373966302,
     -0.04551650597210013,
     -0.04010214944192127,
     -0.03576253848356316,
     -0.0322292075469226,
     -0.029307555889062215,
     -0.026857181875325212
]

lower_bounds = Float64[]
upper_bounds = Float64[]

J1 = 1.0
for (i, J2) in enumerate(0.1:0.2:2.0)
     # i = 1
     # J1, J2 = 1.0, 0.1 + (i - 1) * 0.2

     @ncpolyvar x[1:N] y[1:N] z[1:N]

     obj = one(T) * z[1] * z[2]

     ham = sum(T(J1 / 4) * z[i] * z[mod1(i + 1, N)] + T(J2 / 2) * x[i] for i in 1:N)

     eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

     ineq_cons = [ham - energy_lower_bounds[i] * N]

     pop_l = cpolyopt(obj; eq_constraints=eq_cons, ineq_constraints=ineq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

     pop_u = cpolyopt(-obj; eq_constraints=eq_cons, ineq_constraints=ineq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

     solver_config = SolverConfig(optimizer=SOLVER, order=2)

     single_ineq_cons = [energy_upper_bounds[i] * N - ham]

     res_l = cs_nctssos_with_entry(pop_l, solver_config, single_ineq_cons; dualize=true)
     res_u = cs_nctssos_with_entry(pop_u, solver_config, single_ineq_cons; dualize=true)

     # is_solved_and_feasible(res_l.model) || error("Lower bound problem not solved or feasible at J2 = $J2")
     # is_solved_and_feasible(res_u.model) || error("Upper bound problem not solved or feasible at J2 = $J2")

     # res_l.objective / 4
     # s1s2vals[i]
     # -res_u.objective / 4

     push!(lower_bounds, res_l.objective / 4)
     push!(upper_bounds, -res_u.objective / 4)
end

lower_bounds
upper_bounds

upper_bounds .- lower_bounds
s1s2vals

using DelimitedFiles

open("res.txt","w") do io
     writedlm(io, hcat(lower_bounds, upper_bounds))
end
