using Test, NCTSSoS.FastPolynomials

@testset "Arithmetic" begin
    @testset "Variable Multiplication" begin
        @ncpolyvar x y z
        p1 = 1 * x
        @test p1 isa Polynomial{Float64}
        @test p1.coeffs == [1]
        @test p1.monos == [monomial([x], [1])]
    end

    @testset "Variable Addition" begin
        @ncpolyvar x y z

        p1 = x + y
        @test p1 isa Polynomial{Float64}
        @test p1.coeffs == [1.0, 1.0]
        @test p1.monos == [monomial([x], [1]), monomial([y], [1])]

        p2 = p1 + z
        @test p2.coeffs == [1.0, 1.0, 1.0]
        @test p2.monos == [monomial([x], [1]), monomial([y], [1]), monomial([z], [1])]
    end
    @testset "Monomial multiplication" begin
        @ncpolyvar x y z

        mono1 = monomial([x, y], [1, 2])
        mono2 = monomial([x, z], [3, 4])

        @test mono1 * mono2 == monomial([x, y, x, z], [1, 2, 3, 4])

        mono3 = monomial(Variable[], Int64[])
        @test mono1 * mono3 == mono1
        @test mono3 * mono1 == mono1

        mono4 = monomial([x, y], [1, 0])
        @test mono4 * mono1 == monomial([x, y], [2, 2])
    end

    @testset "Polynomial Addition" begin
        @ncpolyvar x y z
        p1 = Polynomial([1.0, 2.0], [monomial([x], [1]), monomial([y], [2])])
        p_sum1 = p1 + monomial([z], [1])
        @test p_sum1 == Polynomial(
            [1.0, 2.0, 1.0], [monomial([x], [1]), monomial([y], [2]), monomial([z], [1])]
        )
        p2 = monomial([x], [1]) + monomial([y], [2])
        @test p2.coeffs == [1.0, 1.0]
        @test p2.monos == [monomial([x], [1]), monomial([y], [2])]

        p3 = Polynomial(Int[2], [monomial([x, y], [1, 2])])
        p4 = Polynomial(Float32[1.0], [monomial([z], [3])])
        p34 = p3 + p4
        @test p34 isa Polynomial{Float32}
        @test p34.coeffs == [2.0, 1.0]
        @test p34.monos == [monomial([x, y], [1, 2]), monomial([z], [3])]
    end

    @testset "Scaling Polynomial" begin
        @ncpolyvar x y z
        p1 = Polynomial([1.0, 2.0], [monomial([x], [1]), monomial([y], [2])])
        p2 = Float32(2.0) * p1
        @test p2.coeffs == [2.0, 4.0]
        @test p2.monos == [monomial([x], [1]), monomial([y], [2])]
    end

    @testset "Subtraction" begin
        @ncpolyvar x y z
        @test x - y == Polynomial([-1.0, 1.0], [monomial([y], [1]), monomial([x], [1])])

        @test x^2 - y == Polynomial([-1.0, 1.0], [monomial([y], [1]), monomial([x], [2])])

        @test x^2 - y^2 == Polynomial([-1.0, 1.0], [monomial([y], [2]), monomial([x], [2])])

        p1 = 1.0 * x^2
        p2 = 2.0 * y^2
        @test p1 - p2 == Polynomial([-2.0, 1.0], [monomial([y], [2]), monomial([x], [2])])
    end

    @testset "*" begin
        @ncpolyvar x y z

        p1 = 1.0 * x^2
        @test p1 isa Polynomial{Float64}

        p2 = x^2 + y^3
        @test p2 isa Polynomial{Float64}
        @test p2 == Polynomial([1.0, 1.0], [monomial([x], [2]), monomial([y], [3])])

        m1 = x * y
        @test m1.vars == [x, y]
        @test m1.z == [1, 1]

        m2 = x * y * z
        @test m2.vars == [x, y, z]
        @test m2.z == [1, 1, 1]
    end
end
