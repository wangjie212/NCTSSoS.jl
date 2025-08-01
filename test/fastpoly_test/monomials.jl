using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: neat_dot, star, _neat_dot3

@testset "Monomials" begin
    @testset "Creation" begin
        @ncpolyvar x y z
        mono1 = monomial([x, y, x, z], [30, 2, 0, 1])
        @test mono1.vars == [x, y, z]
        @test mono1.z == [30, 2, 1]

        @test_throws ArgumentError monomial([x, y], [1, 2, 3]) # wrong number of variables and exponents

        mono2 = monomial((x, y, z), (1, 2, 3))
        @test mono2.vars == [x, y, z]
        @test mono2.z == [1, 2, 3]

        mono3 = monomial([x, y, z, y, z], [1, 1, 0, 1, 3])
        @test mono3.vars == [x, y, z]
        @test mono3.z == [1, 2, 3]

        @test one(Monomial) == monomial([], [])

        mono4 = monomial(mono3)
        @test mono4 === mono3

        mono5 = monomial(x)
        @test mono5.vars == [x]
    end

    @testset "Utils" begin
        @ncpolyvar x y z

        mono1 = monomial([x, y, z, x], [30, 2, 0, 1])
        mono2 = monomial([x, y, x, z], [0, 0, 0, 0])

        @test degree(mono1) == 33
        @test degree(mono2) == 0

        @test variables(mono1) == [x, y]
        @test variables(mono2) == typeof(x)[]
    end

    @testset "Hash" begin
        @ncpolyvar x y z
        mono1 = monomial([x, y, z], [1, 2, 3])
        mono2 = monomial([x, y, z], [1, 3, 3])
        mono3 = monomial([x, z, y, z], [1, 0, 2, 3])

        @test hash(mono1) != hash(mono2)
        @test hash(mono1) == hash(mono3)

        mono4 = monomial([x], [0])
        @test hash(mono4) == hash(monomial([], []))
    end

    @testset "Star Operation" begin
        @ncpolyvar x y z
        mono1 = monomial([x, y, z], [2, 0, 1])

        # NOTE: I am assuming all variables are Hermitian
        mono1_star = star(mono1)

        @test mono1_star.vars == [z, x]
        @test mono1_star.z == [1, 2]

        mono2 = monomial([x, y, z], [0, 0, 0])
        mono2_star = star(mono2)

        @test isempty(mono2_star.vars)
        @test isempty(mono2_star.z)

        mono3 = monomial([x, y, z], [1, 1, 1])
        mono3_star = star(mono3)
        @test mono3_star.vars == [z, y, x]
        @test mono3_star.z == [1, 1, 1]
    end

    @testset "neat_dot" begin
        @ncpolyvar x y z
        mono1 = monomial([x, y], [1, 0])

        mono2 = monomial([x, y], [1, 1])

        @test neat_dot(mono1, mono2) == monomial([x, y], [2, 1])
        @test neat_dot(mono2, mono2) == monomial([y, x, y], [1, 2, 1])
    end

    @testset "_neat_dot3" begin
        @ncpolyvar x y z
        mono1 = monomial([x, y], [1, 0])
        mono2 = monomial([x, y], [1, 1])
        mono3 = monomial([z, x], [2, 1])

        @test _neat_dot3(mono1, mono2, mono3) == neat_dot(mono1, mono2 * mono3)

        mono1 = monomial([x, y, x], [1, 0, 3])
        mono2 = monomial([x], [2])
        mono3 = monomial([x, z], [2, 1])

        @test _neat_dot3(mono1, mono2, mono3) == neat_dot(mono1, mono2 * mono3)


        mono1 = monomial([x, y], [1, 0])
        mono2 = monomial([], [])
        mono3 = monomial([z, x], [2, 1])

        @test _neat_dot3(mono1, mono2, mono3) == neat_dot(mono1, mono2 * mono3)
    end

    @testset "degree" begin
        @ncpolyvar x y z
        mono1 = monomial([x, y, z], [1, 2, 3])
        mono2 = monomial([x, y], [0, 0])
        mono3 = monomial([x, y, z], [0, 0, 0])

        @test degree(mono1) == 6
        @test degree(mono2) == 0
        @test degree(mono3) == 0
    end
end
