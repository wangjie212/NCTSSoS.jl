using Test, NCTSSoS
using MosekTools

N = 10
@ncpolyvar x[1:N] y[1:N] z[1:N]

ham = sum(1 / 4 * op[i] * op[mod1(i + 1, N)] for op in [x, y, z] for i in 1:N)

pop = PolyOpt(ham; comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=Mosek.Optimizer)

# how do I deal with complex constraint?
cs_nctssos(pop, solver_config)