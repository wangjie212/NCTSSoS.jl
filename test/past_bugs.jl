using Test, NCTSSoS
using MosekTools

@testset "Pauli X Z" begin
    @ncpolyvar x[1:2] z[1:2]

    algebra_constraints = [x[1] * z[1] + z[1] * x[1], x[2] * z[2] + z[2] * x[2]]
    Hamiltonian_objective = 0.5 * x[1] + 0.5 * x[2]
    PXP_constraints = [1.0 * 1 + -1.0 * z[1] + -1.0 * z[2] + 1.0 * z[1]z[2]]
    comm_gps = [[x[i], z[i]] for i in 1:2]

    pop = polyopt(Hamiltonian_objective; eq_constraints=vcat(PXP_constraints, algebra_constraints), comm_gps=comm_gps, is_unipotent=true)

    order = 2
    res = cs_nctssos(pop, SolverConfig(optimizer=Mosek.Optimizer, cs_algo=NoElimination(), ts_algo=NoElimination(), order=order))

    @test res.objective ≈ -sqrt(2) / 2.0
end

@testset "Pauli X Y Z" begin
    N = 2
    @ncpolyvar x[1:2] y[1:2] z[1:2]

    T = ComplexF64

    algebra_constraints = [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N]
    Hamiltonian_objective = zero(T) + 0.5 * x[1] + 0.5 * x[2]
    PXP_constraints = [zero(T) + 1.0 * 1 + -1.0 * z[1] + -1.0 * z[2] + 1.0 * z[1]z[2]]
    comm_gps = [[x[i], y[i], z[i]] for i in 1:2]


    pop = polyopt(Hamiltonian_objective; eq_constraints=[PXP_constraints; algebra_constraints...], comm_gps=comm_gps, is_unipotent=true)

    order = 2
    res = cs_nctssos(pop, SolverConfig(optimizer=Mosek.Optimizer, cs_algo=NoElimination(), ts_algo=NoElimination(), order=order))

    @test res.objective ≈ -sqrt(2) / 2.0
end


@testset "Bug Fix for Issue #125" begin
    @ncpolyvar x[1:2] z[1:2]

    algebra_constraints = [x[1] * z[1] + z[1] * x[1], x[2] * z[2] + z[2] * x[2]]
    Hamiltonian_objective = -1.0 * 1 + 0.5 * z[1] + 0.5 * z[2]
    PXP_constraints = [1.0 * 1 + -1.0 * z[1] + -1.0 * z[2] + 1.0 * z[1]z[2]]
    comm_gps = [[x[i], z[i]] for i in 1:2]

    comm_constraints = [x[1] * x[2] - x[2] * x[1], z[1] * z[2] - z[2] * z[1]]

    pop = polyopt(Hamiltonian_objective; eq_constraints=vcat(PXP_constraints, algebra_constraints), comm_gps=comm_gps, is_unipotent=true)

    order = 4
    res = cs_nctssos(pop, SolverConfig(optimizer=Mosek.Optimizer, cs_algo=NoElimination(), ts_algo=NoElimination(), order=order);dualize=false)


    pop = polyopt(Hamiltonian_objective; eq_constraints=vcat(PXP_constraints, algebra_constraints, comm_constraints), is_unipotent=true)

    order = 2
    res = cs_nctssos(pop, SolverConfig(optimizer=Mosek.Optimizer, cs_algo=NoElimination(), ts_algo=NoElimination(), order=order);dualize=false)
    @test res.objective ≈ -1.
end