using Test, NCTSSoS
using MosekTools

N = 10
@ncpolyvar x[1:N] y[1:N] z[1:N]

ham = sum(ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [x, y, z] for i in 1:N)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

pop = PolyOpt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=Mosek.Optimizer)

# how do I deal with complex constraint?
cs_nctssos(pop, solver_config)

using JuMP

model = GenericModel{Float64}()

@variable(model, y[1:2])

@constraint(model, im * y[1] = y[2])

