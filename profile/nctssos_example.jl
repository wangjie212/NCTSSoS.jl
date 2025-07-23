using BenchmarkTools
using MosekTools, NCTSSoS
using Profile 
using JET

using NCTSSoS.FastPolynomials: monomial, monomials, neat_dot, star, SimplifyAlgorithm, simplify!

@ncpolyvar x[1:6]

@benchmark simplify!(a, sa_proj) setup = (
    a = monomial(x[[6,3,6,1,2,5,4]], [1,1,1,2,1,3,2]);
    sa_proj = SimplifyAlgorithm(comm_gps=[x[1:3], x[4:6]], is_unipotent=false, is_projective=true)
)
# 41 ns

@benchmark simplify!(a, sa_proj) setup = (
    a = monomial(x[[6,3,6,1,2,5,4]], [1,1,1,2,1,3,2]);
    sa_proj = SimplifyAlgorithm(comm_gps=[x[1:3], x[4:6]], is_unipotent=true, is_projective=false)
)
# 20 ns

@benchmark simplify!(a, sa_proj) setup = (
    a = monomial(x[[6,3,6,1,2,5,4]], [1,1,1,2,1,3,2]);
    sa_proj = SimplifyAlgorithm(comm_gps=[x[1:3], x[4:6]], is_unipotent=false, is_projective=false)
)
#40

order = 3
n = 30 
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

br = @benchmark result = cs_nctssos($pop, $solver_config) 

result = cs_nctssos(pop, solver_config)
Profile.clear()
@profile result = cs_nctssos(pop, solver_config)
Profile.print(mincount=100)


# constrain moment matrix
using NCTSSoS
using NCTSSoS: iterate_term_sparse_supp, monomials, get_basis, SimplifyAlgorithm, term_sparsities

@ncpolyvar y[1:10]

sa = SimplifyAlgorithm(comm_gps=[sort(y)], is_unipotent=false, is_projective=false)

basis = get_basis(y, 3)
localizing_basis = [get_basis(y,2), get_basis(y,1)]
init_act_supp = basis[[1,300,20,34,48,59,60]] 
mono = 1.0*y[1]*y[2] + 2.0 * y[2]^3
monos = [mono,mono]

@benchmark term_sparsities($init_act_supp, $monos, $basis, $localizing_basis, MMD(), $sa)

Profile.clear()
@profile for _ in 1:30 term_sparsities(init_act_supp, monos, basis, localizing_basis, MMD(), sa) end
Profile.print(mincount=300)

term_sparsities(init_act_supp, monos, basis, localizing_basis, MMD(), sa)


@benchmark iterate_term_sparse_supp($init_act_supp, $mono, $basis, MMD(), $sa)
iterate_term_sparse_supp(init_act_supp, mono, basis, MMD(), sa)

Profile.clear()
@profile for _ in 1:10 iterate_term_sparse_supp(init_act_supp, mono, basis, MMD(), sa) end
Profile.print(mincount=300)

a = monomial(x[[6, 3, 6, 1, 2, 5, 4]], [1, 1, 1, 2, 1, 3, 2]);
comm_gps = [x[1:3], x[4:6]]
