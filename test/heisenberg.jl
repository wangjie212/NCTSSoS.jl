using Test, NCTSSoS
using MosekTools

SOLVER = Mosek.Optimizer

N = 6
J1 = 1.0
@ncpolyvar x[1:N] y[1:N] z[1:N]

ham = sum(ComplexF64(J1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [x, y, z] for i in 1:N)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

pop = PolyOpt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=SOLVER, mom_order=2, ts_algo = MMD())

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(pop,res,solver_config)
@test res.objective / N ≈ -0.467129 atol=1e-6


N = 6 
J1 = 1.0
J2 = 0.2
@ncpolyvar x[1:N] y[1:N] z[1:N]

ham = sum(ComplexF64(J1 / 4) * op[i] * op[mod1(i + 1, N)] + ComplexF64(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

pop = PolyOpt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=SOLVER, mom_order=2, ts_algo = MMD())

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(pop,res,solver_config)

@test res.objective / N ≈ -0.4270083225302217 atol=1e-6


N = 6
J1 = 1.0
energy_upper_bounds = [-0.46712927295533224, -0.4568115059758944, -0.44667162694019047, -0.43672967371831345, -0.42700832253022164, -0.41753322039150603, -0.40833333333333327, -0.3994412953910567, -0.3908937341178951, -0.3827315360410695, -0.37499999999999983]
s0s1_vals = [-0.1557097576496684, -0.15568115487708442, -0.1555868522576389, -0.15541222562212834, -0.15513979550197515, -0.15474876760316292, -0.1542145593761434, -0.153508352433215, -0.15259673575852825, -0.15144153948520134, -0.16919727528839815]
for (idx,J2) in enumerate(0.0:0.05:0.5)
@ncpolyvar x[1:N] y[1:N] z[1:N]

ham = sum(ComplexF64(J1 / 4) * op[i] * op[mod1(i + 1, N)] + ComplexF64(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

pop = PolyOpt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=SOLVER, mom_order=2, ts_algo = MMD())

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(pop,res,solver_config)

energy_lower_bound = res.objective 

@ncpolyvar x[1:N] y[1:N] z[1:N]

obj = one(ComplexF64) * x[1] * x[2]

ham = sum(ComplexF64(J1 / 4) * op[i] * op[mod1(i + 1, N)] + ComplexF64(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

ineq_cons = [energy_upper_bounds[idx] - ham, ham - energy_lower_bound]

pop = PolyOpt(obj; eq_constraints=eq_cons, ineq_constraints=ineq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=SOLVER, mom_order=2, ts_algo = MMD())

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(pop,res,solver_config)

operator_lower_bound = res.objective

pop = PolyOpt(-obj; eq_constraints=eq_cons, ineq_constraints=ineq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=SOLVER, mom_order=2, ts_algo = MMD())

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(pop,res,solver_config)

operator_upper_bound = -res.objective
end

Nx = 3
Ny = 3
N = Nx * Ny
J1 = 1.0
J2 = 0.0
@ncpolyvar x[1:N] y[1:N] z[1:N]

LI = LinearIndices((1:Nx, 1:Ny))
LI[CartesianIndex(1, 2)]

ham = sum(ComplexF64(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(i, mod1(j + 1, Ny))]] + ComplexF64(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), j)]] + ComplexF64(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j + 1, Ny))]] + ComplexF64(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j - 1, Ny))]] for op in [x, y, z] for i in 1:Nx for j in 1:Ny)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

pop = PolyOpt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=SOLVER, mom_order=3, cs_algo=MF(), ts_algo=MMD())

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(pop,res,solver_config) 
res = cs_nctssos_higher(pop,res,solver_config) 


res = cs_nctssos_higher(pop,res,solver_config) # -4.390300714054776 = -0.4878111904505307 * N
res = cs_nctssos_higher(pop,res,solver_config) # -4.381164563801521 = -0.48679606264461345 

@test res.objective / N ≈ -0.44100019443650207 atol = 1e-6



# DMRG Reference Code
# using MPSKitModels
# using MPSKit
# using KrylovKit
# using TensorKit

# N = 6
# χ = 100
# max_bond_dimension = ℂ^χ
# physical_space = ℂ^2
# state = FiniteMPS(rand, ComplexF64, N, physical_space, max_bond_dimension)
# lattice = FiniteChain(N)

# h = 0.25
# h2 = 0.0
# # h2 = 0.25*0.2
# H_heisenberg = @mpoham sum(h * op(){lattice[i],lattice[mod1(i + 1, N)]} + h2 * op(){lattice[i],lattice[mod1(i + 2, N)]} for op in (S_xx, S_yy, S_zz) for i in 1:N)

# ground_state, cache, delta = find_groundstate(state, H_heisenberg, DMRG())

# real(expectation_value(ground_state, H_heisenberg)) / N * 4

# N = 6
# χ = 100
# max_bond_dimension = ℂ^χ
# physical_space = ℂ^2
# state = FiniteMPS(rand, ComplexF64, N, physical_space, max_bond_dimension)
# lattice = FiniteChain(N)

# h = 0.25
# h2 = 0.25*0.2
# H_heisenberg = @mpoham sum(h * op(){lattice[i],lattice[mod1(i + 1, N)]} + h2 * op(){lattice[i],lattice[mod1(i + 2, N)]} for op in (S_xx, S_yy, S_zz) for i in 1:N)

# ground_state, cache, delta = find_groundstate(state, H_heisenberg, DMRG())

# real(expectation_value(ground_state, H_heisenberg)) / N * 4


# Nx = 3
# Ny = 3
# N = Nx * Ny 
# χ = 500
# max_bond_dimension = ℂ^χ
# physical_space = ℂ^2
# state = FiniteMPS(rand, ComplexF64, N, physical_space, max_bond_dimension)
# lattice = FiniteCylinder(Nx, N)


# h = 0.25 
# h2 = 0.0
# H_heisenberg = @mpoham sum(h * op(){lattice[i, j],lattice[mod1(i + 1, Nx), j]} + h * op(){lattice[i, j],lattice[i, mod1(j + 1, Ny)]} + h2 * op(){lattice[i, j],lattice[mod1(i + 1, Nx), mod1(j + 1, Ny)]} + h2 * op(){lattice[i, j],lattice[mod1(i + 1, Nx), mod1(j - 1, Ny)]} for op in (S_xx, S_yy, S_zz) for i in 1:Nx for j in 1:Ny)

# ground_state, cache, delta = find_groundstate(state, H_heisenberg, DMRG())
# # -0.44100019443650207
# real(expectation_value(ground_state, H_heisenberg)) / N * 4
# # at Nx , Ny = 4,4 result matched wiht Certifying Ground State Table IX



