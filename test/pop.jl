using Test, NCTSSoS, NCTSSoS.FastPolynomials

@testset "PolyOpt Constructor" begin
    nvars = 10
    ncons = 3
    @ncpolyvar x[1:nvars]
    objective = 1.0 * sum(x .^ 2)
    constraints = [1.0 * sum(i .* x) for i = 1:ncons]

    @testset "Unconstrained" begin
        pop = PolyOpt(objective)

        @test isempty(pop.eq_constraints)
        @test isempty(pop.ineq_constraints)

        @test sort(pop.variables) == sort(x)
        @test sort.(pop.comm_gps) == sort.([x])
        @test !pop.is_unipotent
        @test !pop.is_projective
    end

    @testset "Constrainted Optimization Problem" begin
        pop = PolyOpt(objective; ineq_constraints = constraints)

        @test pop.ineq_constraints == constraints
        @test isempty(pop.eq_constraints)

        pop = PolyOpt(objective; ineq_constraints = [constraints; sum(x)])

        @test length(pop.ineq_constraints) == ncons

        pop = PolyOpt(
            objective;
            eq_constraints = constraints[2:2:end],
            ineq_constraints = constraints[1:2:end],
            is_unipotent = false,
            is_projective = true,
        )

        @test length(pop.eq_constraints) == 1
        @test length(pop.ineq_constraints) == 2
        @test pop.is_unipotent == false
        @test pop.is_projective == true
    end

    @testset "Invalid Input" begin
        @test_throws AssertionError PolyOpt(
            objective;
            is_unipotent = true,
            is_projective = true,
        )
        p1 = Polynomial(
            [1, 1],
            [monomial([x[1], x[2]], [1, 1]), monomial([x[2], x[3]], [1, 1])],
        )
        @test_throws AssertionError PolyOpt(p1)
        @ncpolyvar y[1:nvars]
        @test_throws AssertionError PolyOpt(objective; comm_gps = [[x], [y]])
    end
end


@testset "StatePolyOpt Constructor" begin
    @testset "Example 7.2.1" begin
        @ncpolyvar x[1:2] y[1:2]
        sp1 = sum([1.0, 1.0] .* map(a -> prod(ς.(a)), [[x[1] * y[2]], [y[2] * x[1]]]))
        sp2 = sum([1.0, -1.0] .* map(a -> prod(ς.(a)), [[x[1] * y[1]], [x[2] * y[2]]]))
        sp = sum([sp1 * sp1, sp2 * sp2])

        sp1_sq = sum(
            [1.0, 1.0, 1.0, 1.0] .* map(
                a -> prod(ς.(a)),
                [
                    [x[1] * y[2], x[1] * y[2]],
                    [y[2] * x[1], y[2] * x[1]],
                    [x[1] * y[2], y[2] * x[1]],
                    [y[2] * x[1], x[1] * y[2]],
                ],
            ),
        )
        sp2_sq = sum(
            [1.0, -1.0, -1.0, 1.0] .* map(
                a -> prod(ς.(a)),
                [
                    [x[1] * y[1], x[1] * y[1]],
                    [x[1] * y[1], x[2] * y[2]],
                    [x[2] * y[2], x[1] * y[1]],
                    [x[2] * y[2], x[2] * y[2]],
                ],
            ),
        )
        true_obj = sum([sp1_sq, sp2_sq])

        pop = PolyOpt(sp * one(Monomial); is_unipotent = true, comm_gps = [x, y])
        @test pop.objective == true_obj * one(Monomial)
        @test isempty(pop.eq_constraints)
        @test isempty(pop.ineq_constraints)
        @test pop.is_unipotent == true
        @test pop.is_projective == false
        @test pop.comm_gps == [x, y]
    end
    @testset "Example 7.2.2" begin
        @ncpolyvar x[1:3] y[1:3]
        cov(a, b) = 1.0 * ς(x[a] * y[b]) - ς(x[a]) * ς(y[b])
        sp =
            cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) +
            cov(3, 1) - cov(3, 2)

        pop = PolyOpt(sp * one(Monomial); is_unipotent = true, comm_gps = [x, y])
        true_obj = sum(
            [
                1.0,
                -1.0,
                1.0,
                -1.0,
                1.0,
                -1.0,
                1.0,
                -1.0,
                1.0,
                -1.0,
                -1.0,
                1.0,
                1.0,
                -1.0,
                -1.0,
                1.0,
            ] .* map(
                a -> prod(ς.(a)),
                ([
                    [x[1] * y[1]],
                    [x[1], y[1]],
                    [x[1] * y[2]],
                    [x[1], y[2]],
                    [x[1] * y[3]],
                    [x[1], y[3]],
                    [x[2] * y[1]],
                    [x[2], y[1]],
                    [x[2] * y[2]],
                    [x[2], y[2]],
                    [x[2] * y[3]],
                    [x[2], y[3]],
                    [x[3] * y[1]],
                    [x[3], y[1]],
                    [x[3] * y[2]],
                    [x[3], y[2]],
                ]),
            ),
        )
        @test pop.objective == true_obj * one(Monomial)
        @test isempty(pop.eq_constraints)
        @test isempty(pop.ineq_constraints)
        @test pop.is_unipotent == true
        @test pop.comm_gps == [x, y]

        @ncpolyvar z[1:3]
        @test_throws AssertionError PolyOpt(sp*one(Monomial); comm_gps = [x, y, z])
    end

    @testset "Example 8.1.2" begin
        @ncpolyvar A[1:3] B[1:3]
        J1 = 0.5 * (ς(A[1]) + ς(A[2]) + ς(A[3]) + ς(B[1] * A[1]) + ς(B[1] * A[2]))

        J2 =
            0.5 * (ς(A[1]) + ς(A[2]) - ς(A[3]) + ς(B[2] * A[1]) - ς(B[2] * A[2])) +
            0.5 * (
                ς(A[1] * B[3] * A[1]) - ς(A[1] * B[3] * A[2]) - ς(A[2] * B[3] * A[1]) +
                ς(A[2] * B[3] * A[2])
            )


        L = 4.0 + ς(A[1]) + ς(A[2])

        sp = sum([
            2.0 * J1 * J2,
            2.0 * J1 * L,
            2.0 * J2 * L,
            -1.0 * J1 * J1,
            -1.0 * J2 * J2,
            -1.0 * L * L,
        ])
        pop = PolyOpt(sp*one(Monomial); is_unipotent = true, comm_gps = [A, B])
        @test isempty(pop.eq_constraints)
        @test isempty(pop.ineq_constraints)
        @test pop.is_unipotent == true
        @test pop.comm_gps == [A, B]
    end
end
