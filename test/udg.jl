using Graphs, NCTSSoS
using DynamicPolynomials
using MosekTools
using JuMP

n = 18
L = 5.0
cutoff= 1.0
order = 2

# 2d case with euclidean norm
# g, edgs = euclidean_graph(n, 2; seed=110, L=L, p=2.0, cutoff=cutoff, bc=:open)

g = star_graph(n) 

@polyvar x[1:n] y[1:n] 

ijpairs = [(i, j) for i in 1:n for j in 1:n if j > i]

ijpair2idx = Dict(zip(ijpairs, collect(1:length(ijpairs))))


f = 1.0* polynomial(one(x[1])) 
# f = sum(vec([(ispos[ijpair2idx[(i, j)]] - 1) * (has_edge(g, i, j) ? (cutoff^2 - (x[i] - x[j])^2 - (y[i] - y[j])^2) : ((x[i] - x[j])^2 + (y[i] - y[j])^2) - cutoff^2) for i in 1:n for j in 1:n if j > i]))

cons = [
    [
        (cutoff^2 - (x[edg.src] - x[edg.dst])^2 - (y[edg.src] - y[edg.dst])^2)
        for edg in edges(g)
    ];[
        ((x[i] - x[j])^2 + (y[i] - y[j])^2 - cutoff^2)
        for i in 1:n, j in 1:n if j > i && !(has_edge(g, i, j))
    ]
]

pop = PolyOpt(f; constraints=cons, is_equality=[false for _ in cons])

cs_algo = NoElimination() 
ts_algo = MMD()

using NCTSSoS: correlative_sparsity, iterate_term_sparse_supp, sorted_union, symmetric_canonicalize, neat_dot, moment_relax, sos_dualize
corr_sparsity = correlative_sparsity(pop, order, cs_algo)

cliques_objective = [reduce(+, [issubset(effective_variables(mono), clique) ? coef * mono : zero(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

# prepare the support for each term sparse localizing moment
initial_activated_supp = [
    sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, pop.constraints[cons_idx];init=[one(pop.variables[1])]), [neat_dot(b, b) for b in idcs_bases[1]])
    for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)
]

cliques_term_sparsities = map(zip(initial_activated_supp, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
    [iterate_term_sparse_supp(activated_supp, poly, basis, ts_algo) for (poly, basis) in zip([one(pop.objective); pop.constraints[cons_idx]], idcs_bases)]
end

moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, corr_sparsity.global_cons, cliques_term_sparsities)
# sos_problem = sos_dualize(moment_problem)

set_optimizer(moment_problem.model, Mosek.Optimizer)
# set_optimizer(sos_problem.model, Mosek.Optimizer)

# optimize!(sos_problem.model)
optimize!(moment_problem.model)
# is_solved_and_feasible(sos_problem.model)
# solution_summary(sos_problem.model)
is_solved_and_feasible(moment_problem.model)
solution_summary(moment_problem.model)


objective_value(sos_problem.model)

# star_graph(7) -2.43 e -13


ys_idx = [2, 3, 4, 5, 6, 7, 8]
xs_idx = [9, 10, 11, 12, 13, 14, 15]

xs = value.(moment_problem.model[:y])[xs_idx]
ys = value.(moment_problem.model[:y])[ys_idx]

for (i,a) in enumerate(value.(moment_problem.model[:y]))
    if !isapprox(a, 0.0, atol=1e-6)
        @show i a
    end
end
a = 1
b = 2

# problem is (<a>-<b>)^2 is not equal to <a^2> - 2<ab> + <b^2>
for i in 1:n, j in (i+1):n
    xsq = value(moment_problem.monomap[x[i]^2]) - 2 * value(moment_problem.monomap[x[i]*x[j]]) + value(moment_problem.monomap[x[j]^2])
    xsq2 = (value(moment_problem.monomap[x[i]])-value(moment_problem.monomap[x[j]]))^2
    ysq = value(moment_problem.monomap[y[i]^2]) - 2 * value(moment_problem.monomap[y[i]*y[j]]) + value(moment_problem.monomap[y[j]^2])
    ysq2 = (value(moment_problem.monomap[y[i]])-value(moment_problem.monomap[y[j]]))^2
    @assert isapprox(xsq, xsq2, atol=1e-6)
    @assert isapprox(ysq, ysq2, atol=1e-6)
    if Edge(i, j) in edges(g_not)
        @assert xsq + ysq <= cutoff^2
    else
        @assert xsq + ysq >= cutoff^2
    end
end
