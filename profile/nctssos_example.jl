using BenchmarkTools
using MosekTools, NCTSSoS
using Profile 
using JET

order = 3
n = 20 
@ncpolyvar x[1:n]
f = 0.0
for i = 1:n
    jset = max(1, i - 5):min(n, i + 1)
    jset = setdiff(jset, i)
    f += (2x[i] + 5 * x[i]^3 + 1)^2
    f -= sum([
        4x[i] * x[j] +
        10x[i]^3 * x[j] +
        2x[j] +
        4x[i] * x[j]^2 +
        10x[i]^3 * x[j]^2 +
        2x[j]^2 for j in jset
    ])
    f += sum([
        x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset
    ])
end

cons = vcat([(1 - x[i]^2) for i = 1:n], [(x[i] - 1 / 3) for i = 1:n])

pop = polyopt(f; ineq_constraints=cons)
solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=order,
    cs_algo=MF(), ts_algo=MMD())


# @benchmark result = cs_nctssos($pop, $solver_config) 

result = cs_nctssos(pop, solver_config)
Profile.clear()
@profile result = cs_nctssos(pop, solver_config)
Profile.print(mincount=300)


# constrain moment matrix
using NCTSSoS
using NCTSSoS: iterate_term_sparse_supp, monomials, get_basis, SimplifyAlgorithm

@ncpolyvar y[1:10]

sa = SimplifyAlgorithm(comm_gps=[y], is_unipotent=false, is_projective=false)

basis = get_basis(y, 3)
init_act_supp = basis[[1,300,20,34,48,59,60]] 
mono = 1.0*y[1]*y[2] + 2.0 * y[2]^3

@benchmark iterate_term_sparse_supp($init_act_supp, $mono, $basis, MMD(), $sa)

Profile.clear()
@profile for _ in 1:10 iterate_term_sparse_supp(init_act_supp, mono, basis, MMD(), sa) end
Profile.print(mincount=300)

@code_warntype iterate_term_sparse_supp(init_act_supp, mono, basis, MMD(), sa)
report = @report_opt iterate_term_sparse_supp(init_act_supp, mono, basis, MMD(), sa)
report = @report_call iterate_term_sparse_supp(init_act_supp, mono, basis, MMD(), sa)

show(report)

report

ts = iterate_term_sparse_supp(init_act_supp, mono, basis, MMD(), sa)

using NCTSSoS.FastPolynomials: neat_dot, star , monomial, Monomial, ProdMonomial
using NCTSSoS

using BenchmarkTools
using MosekTools, NCTSSoS
using Profile 
using JET

@ncpolyvar x[1:10]

@benchmark Monomial(vars, expos) setup = (vars = x[[[6, 8, 7, 1, 2, 5, 8, 3, 4, 2]; [5, 1, 3, 7, 4, 8, 7, 6, 3, 9]]]; expos = [8, 5, 3, 6, 5, 10, 2, 8, 10, 7, 8, 4, 3, 2, 9, 8, 10, 6, 8, 9])


@benchmark star(a) setup = (a = monomial(x[[6, 8, 7, 1, 2, 5, 8, 3, 4, 2]], [8, 5, 3, 6, 5, 10, 2, 8, 10, 7]))

@benchmark a * b setup = (a = monomial(x[[6, 8, 7, 1, 2, 5, 8, 3, 4, 2]], [8, 5, 3, 6, 5, 10, 2, 8, 10, 7]); b = monomial(x[[5, 1, 3, 7, 4, 8, 7, 6, 3, 9]], [8, 4, 3, 2, 9, 8, 10, 6, 8, 9]))

@benchmark neat_dot(a, b) setup = (a = monomial(x[[6, 8, 7, 1, 2, 5, 8, 3, 4, 2]], [8, 5, 3, 6, 5, 10, 2, 8, 10, 7]); b = monomial(x[[5, 1, 3, 7, 4, 8, 7, 6, 3, 9]], [8, 4, 3, 2, 9, 8, 10, 6, 8, 9]))

using NCTSSoS.FastPolynomials: _neat_dot3
@benchmark neat_dot(a, b*c) setup=(a = monomial(x[[6, 8, 7, 1, 2, 5, 8, 3, 4, 2]], [8, 5, 3, 6, 5, 10, 2, 8, 10, 7]); b = monomial(x[[5, 1, 3, 7, 4, 8, 7, 6, 3, 9]], [8, 4, 3, 2, 9, 8, 10, 6, 8, 9]); c = monomial(x[[3, 4, 5, 7]],[4, 4, 2, 4]))

@benchmark _neat_dot3(a,b,c) setup=(a = monomial(x[[6, 8, 7, 1, 2, 5, 8, 3, 4, 2]], [8, 5, 3, 6, 5, 10, 2, 8, 10, 7]); b = monomial(x[[5, 1, 3, 7, 4, 8, 7, 6, 3, 9]], [8, 4, 3, 2, 9, 8, 10, 6, 8, 9]); c = monomial(x[[3, 4, 5, 7]],[4, 4, 2, 4]))


# neat_dot (168ns) = product (75ns) + reverse etc  + monomial creation (43ns)