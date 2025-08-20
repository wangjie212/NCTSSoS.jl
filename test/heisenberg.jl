using NCTSSoS, NCTSSoS.FastPolynomials, Test
using MosekTools
using JuMP

SOLVER = Mosek.Optimizer

@testset "Bounding s0s1" begin
    T = ComplexF64
    N = 4
    J1 = 1.0
    energy_upper_bounds = [
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

    energy_lower_bounds = [
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

    op_lower_bounds = Float64[]
    op_upper_bounds = Float64[]
    # for (idx, J2) in enumerate(0.0:0.05:2.0)
        idx, J2 = 1 , 0.1
        @ncpolyvar x[1:N] y[1:N] z[1:N]

        obj = one(T) * x[1] * x[2]

        ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] + T(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

        eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

        push!(eq_cons, obj - 0.5)

        # because energy bounds are negative
        # ineq_cons = [obj-0.5]
        # ineq_cons = [0.9 * N * energy_upper_bounds[idx] - ham, ham - energy_lower_bounds[idx] * N * 1.1]

        # why is negative upper bound of hamiltonian not able to be handled?
        # N * energy_lower_bounds[idx]
        # 0.9 * N * energy_upper_bounds[idx]
        # 1.1 * N * energy_lower_bounds[idx]
        # ineq_cons = [0.99 * N * energy_lower_bounds[idx] - ham, ]

        using SCS
        SOLVER = SCS.Optimizer

        pop = cpolyopt(obj; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)
        # pop = cpolyopt(obj; eq_constraints=eq_cons, ineq_constraints=ineq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

        solver_config = SolverConfig(optimizer=SOLVER, order=2)

        res = cs_nctssos(pop, solver_config);


        # why is it terminating due to slow progress  
        termination_status(res.model)
        solution_summary(res.model)

        @assert is_solved_and_feasible(res.model)

        res.objective

        push!(op_lower_bounds,res.objective)

        pop = cpolyopt(-obj; eq_constraints=eq_cons, ineq_constraints=ineq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

        solver_config = SolverConfig(optimizer=SOLVER, order=2)

        res = cs_nctssos(pop, solver_config)

        @assert is_solved_and_feasible(res.model)

        push!(op_upper_bounds,-res.objective)
    # end

    println("lower")

    for lb in op_lower_bounds
        println(lb / 4, ",")
    end

    println("upper")

    for lb in op_upper_bounds
        println(lb / 4, ",")
    end
end


@testset "J1 J2 Model" begin
    T = ComplexF64
	N = 4
    J1 = 1.0
    J2s = 0.1:0.1:2.1
    @ncpolyvar x[1:N] y[1:N] z[1:N]

    energy_lower_bounds = zeros(length(J2s))

    for (idx, J2) in enumerate(J2s)
        ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] + T(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

        eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

        pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

        solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())

        res = cs_nctssos(pop, solver_config)

        res = cs_nctssos_higher(pop, res, solver_config)
        energy_lower_bounds[idx] = res.objective / N
    end
    for val in energy_lower_bounds
        println(val, ",")
    end
end

# s0s1_vals = [
#     -0.16666666666666657,
#     -0.16666666666666657,
#     -0.08004709661437667,
#     -2.7856520089408856e-17,
#     -3.408084308945605e-17,
#     -3.899516274508277e-17,
#     6.385764772817207e-17,
#     9.579194098506272e-17,
#     1.1493569501470919e-16,
#     1.2608731144223405e-16,
#     -1.1373213176545312e-16]


@testset "XXX Model" begin
    T = ComplexF64
    N = 6
    J1 = 1.0
    @ncpolyvar x[1:N] y[1:N] z[1:N]

    ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [x, y, z] for i in 1:N)

    eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

    pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

    solver_config = SolverConfig(optimizer=SOLVER, order=2)

    res = cs_nctssos(pop, solver_config)

    @test res.objective / N ≈ -0.467129 atol = 1e-6
end


@testset "J1 J2 Model" begin
	N = 6
    J1 = 1.0
    J2 = 0.2
    @ncpolyvar x[1:N] y[1:N] z[1:N]

    ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] + T(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

    eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

    pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

    solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())

    res = cs_nctssos(pop, solver_config)

    res = cs_nctssos_higher(pop, res, solver_config)

    @test res.objective / N ≈ -0.4270083225302217 atol = 1e-6
end

@testset "2D Model" begin
    Nx = 3
    Ny = 3
    N = Nx * Ny
    J1 = 1.0
    J2 = 0.0
    @ncpolyvar x[1:N] y[1:N] z[1:N]

    LI = LinearIndices((1:Nx, 1:Ny))

    ham = sum(ComplexF64(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(i, mod1(j + 1, Ny))]] + ComplexF64(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), j)]] + ComplexF64(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j + 1, Ny))]] + ComplexF64(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j - 1, Ny))]] for op in [x, y, z] for i in 1:Nx for j in 1:Ny)

    eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

    pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

    solver_config = SolverConfig(optimizer=SOLVER, order=3, cs_algo=MF(), ts_algo=MMD())

    res = cs_nctssos(pop, solver_config)

    res = cs_nctssos_higher(pop, res, solver_config)
    res = cs_nctssos_higher(pop, res, solver_config)

    res = cs_nctssos_higher(pop, res, solver_config) # -4.390300714054776 = -0.4878111904505307 * N
    res = cs_nctssos_higher(pop, res, solver_config) # -4.381164563801521 = -0.48679606264461345 

    @test res.objective / N ≈ -0.44100019443650207 atol = 1e-6
end

# using Pkg
# Pkg.activate("/home/yushengzhao/LearnTNTheHardWay.jl")
# # DMRG Reference Code
# using MPSKitModels
# using MPSKit
# using KrylovKit
# using TensorKit

# N = 4
# χ = 200
# max_bond_dimension = ℂ^χ
# physical_space = ℂ^2
# state = FiniteMPS(rand, ComplexF64, N, physical_space, max_bond_dimension)
# lattice = FiniteChain(N)

# h = 0.25
# h2s = collect(0.1:0.2:2.1) ./4
# upper_bounds = zeros(length(h2s))
# s0s1_vals = zeros(length(h2s))
# for (i, h2) in enumerate(h2s)
#     H_heisenberg = @mpoham sum(h * op(){lattice[i],lattice[mod1(i + 1, N)]} + h2 * op(){lattice[i],lattice[mod1(i + 2, N)]} for op in (S_xx, S_yy, S_zz) for i in 1:N)

#     ground_state, cache, delta = find_groundstate(state, H_heisenberg, DMRG())

#     s0s1_vals[i] = real(expectation_value(ground_state, @mpoham S_xx(){lattice[1],lattice[2]}))

#     upper_bounds[i] = real(expectation_value(ground_state, H_heisenberg)) / N * 4
# end

# for val in upper_bounds
#     println(val, ",")
# end

# for val in s0s1_vals
#     println(val, ",")
# end

# s0s1_vals

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



