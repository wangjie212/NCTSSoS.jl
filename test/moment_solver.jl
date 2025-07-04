using Test, NCTSSoS
using NCTSSoS.FastPolynomials
using JuMP

if Sys.isapple()
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end
using Graphs

using NCTSSoS.FastPolynomials: get_basis 
using NCTSSoS: substitute_variables

@testset "Special Constraint Type " begin
    @testset "CHSH Inequality" begin
        @ncpolyvar x[1:2]
        @ncpolyvar y[1:2]

        f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
        pop = polyopt(f; comm_gps = [x, y], is_unipotent = true)

        solver_config = SolverConfig(optimizer = SOLVER; order = 1)

        result = cs_nctssos(pop, solver_config; dualize=false)

        @test isapprox(result.objective, -2.8284271321623193, atol = 1e-6)
    end

    @testset "I_3322 inequality" begin
        @ncpolyvar x[1:3]
        @ncpolyvar y[1:3]
        f =
            1.0 * x[1] * (y[1] + y[2] + y[3]) +
            x[2] * (y[1] + y[2] - y[3]) +
            x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]
        pop = polyopt(-f; comm_gps = [x, y], is_projective = true)

        solver_config = SolverConfig(optimizer = SOLVER; order = 3)

        result = cs_nctssos(pop, solver_config)

        @test isapprox(result.objective, -0.2508753049688358, atol = 1e-6)
    end
end

@testset "CS TS Example" begin
    order = 3
    n = 10
    @ncpolyvar x[1:n]
    f = 0.0
    for i = 1:n
        jset = max(1, i-5):min(n, i+1)
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

    # NOTE: should not symmetric canonicalize the polynomial since it might
    # take away some support in the polynomial
    cons = vcat([(1 - x[i]^2) for i = 1:n], [(x[i] - 1 / 3) for i = 1:n])

    pop = polyopt(f; ineq_constraints = cons)

    solver_config =
        SolverConfig(optimizer=SOLVER; order=order, cs_algo=MF(), ts_algo=MMD())

    result = cs_nctssos(pop, solver_config; dualize=false)

    @test isapprox(result.objective, 3.011288, atol = 1e-4)
end

@testset "Moment Method Heisenberg Model on Star Graph" begin
    num_sites = 10
    g = star_graph(num_sites)

    true_ans = -1.0

    vec_idx2ij = [(i, j) for i = 1:num_sites for j = (i+1):num_sites]

    findvaridx(i, j) = findfirst(x -> x == (i, j), vec_idx2ij)

    @ncpolyvar pij[1:length(vec_idx2ij)]

    objective = sum(1.0 * pij[[findvaridx(ee.src, ee.dst) for ee in edges(g)]])

    gs = unique!([
        (
            pij[findvaridx(sort([i, j])...)] * pij[findvaridx(sort([j, k])...)] +
            pij[findvaridx(sort([j, k])...)] * pij[findvaridx(sort([i, j])...)] -
            pij[findvaridx(sort([i, j])...)] - pij[findvaridx(sort([j, k])...)] -
            pij[findvaridx(sort([i, k])...)] + 1.0
        ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
        (i != j && j != k && i != k)
    ])


    pop = polyopt(
        objective;
        eq_constraints = gs,
        is_unipotent = true,
    )
    order = 1
    cs_algo = MF()

    solver_config = SolverConfig(
        optimizer = SOLVER,
        order = order,
        cs_algo = cs_algo
    )

    result = cs_nctssos(pop, solver_config; dualize = false)
    @test isapprox(result.objective, true_ans, atol = 1e-6)
end

@testset "Replace DynamicPolynomials variables with JuMP variables" begin
    @ncpolyvar x y z
    poly = 1.0 * x^2 - 2.0 * x * y - 1

    model = GenericModel{Float64}()
    @variable(model, jm[1:13])

    monomap = Dict(get_basis([x, y, z], 2) .=> jm)

    @test substitute_variables(poly, monomap) == 1.0 * jm[7] - 2.0 * jm[5] - jm[1]
end

@testset "Moment Method Example 1" begin
    order = 2
    n = 3
    @ncpolyvar x[1:n]
    f =
        x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] +
        2.0x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6x[3]^2 +
        9x[2]^2 * x[3] +
        9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    pop = polyopt(f)

    @testset "Dense" begin

        solver_config = SolverConfig(
            optimizer = SOLVER,
            order = order,
            cs_algo = NoElimination(),
        )

        result = cs_nctssos(pop, solver_config; dualize = false)

        # NOTE: differs from original test case value since that one is a relaxed in terms of sparsity
        # This value here is obtained by running the master branch with no sparsity relaxation
        @test isapprox(
            result.objective,
            4.372259295498716e-10,
            atol = 1e-6,
        )
    end

    @testset "Sprase" begin
        solver_config = SolverConfig(
            optimizer = SOLVER,
            order = order,
            ts_algo = MMD(),
        )

        result = cs_nctssos(pop, solver_config; dualize = false)

        @test isapprox(result.objective, -0.0035512, atol = 1e-7)
    end
end

@testset "Moment Method Example 2" begin
    order = 2
    n = 2
    @ncpolyvar x[1:2]
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0
    pop = polyopt(f; eq_constraints=[h1], ineq_constraints=[g])


    @testset "Dense" begin
        solver_config = SolverConfig(
            optimizer = SOLVER,
            order = order)

        result = cs_nctssos(pop, solver_config; dualize = false)
        @test isapprox(result.objective, -1.0, atol = 1e-6)
    end

    @testset "Term Sparse" begin
        solver_config = SolverConfig(
            optimizer = SOLVER,
            order = order,
            cs_algo = MF(),
            ts_algo = MMD(),
        )

        result = cs_nctssos(pop, solver_config; dualize = false)

        @test isapprox(result.objective, -1.0, atol = 1e-6)
    end
end

@testset "Moment Method Correlative Sparsity" begin
    n = 3
    @ncpolyvar x[1:n]
    f =
        x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0x[2]^2 - 2x[1] * x[2] * x[1] +
        2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
        6.0 * x[3]^2 +
        9x[2]^2 * x[3] +
        9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

    cons = vcat([1.0 - x[i]^2 for i = 1:n], [x[i] - 1.0 / 3 for i = 1:n])
    order = 3

    pop = polyopt(f; ineq_constraints = cons)

    @testset "Correlative Sparse" begin
        solver_config = SolverConfig(
            optimizer = SOLVER,
            order = order,
            cs_algo = MF(),
        )
        result = cs_nctssos(pop, solver_config; dualize = false)

        # FIXME: reduced accuracy
        # @test is_solved_and_feasible(moment_problem.model)
        @test isapprox(
            result.objective,
            0.9975306427277915,
            atol = 1e-5,
        )
    end

    @testset "Term Sparse" begin
        solver_config = SolverConfig(
            optimizer = SOLVER,
            order = order,
            ts_algo = MMD(),
        )

        result = cs_nctssos(pop, solver_config; dualize = false)

        result = cs_nctssos_higher(pop, result, solver_config;dualize=false)

        @test isapprox(
            result.objective,
            0.9975306427277915,
            atol = 1e-5,
        )
    end
end
