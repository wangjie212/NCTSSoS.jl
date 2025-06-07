using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: degree, Monomial

@testset "Monomials" begin
    @testset "Creation" begin
        @ncpolyvar x y z
        mono1 = Monomial([x, y, x, z], [30, 2, 0, 1])
        @test mono1.vars == [x, y, z]
        @test mono1.z == [30, 2, 1]

        @test_throws ArgumentError Monomial([x, y], [1, 2, 3]) # wrong number of variables and exponents


        mono2 = Monomial((x, y, z), (1, 2, 3))
        @test mono2.vars == [x, y, z]
        @test mono2.z == [1, 2, 3]
    end

    @testset "Utils" begin
        @ncpolyvar x y z

        mono1 = Monomial([x, y, z, x], [30, 2, 0, 1])
        mono2 = Monomial([x, y, x, z], [0, 0, 0, 0])

        @test degree(mono1) == 33
        @test degree(mono2) == 0

        @test variables(mono1) == [x,y]
        @test variables(mono2) == typeof(x)[]
    end

    @testset "Hash" begin
        @ncpolyvar x y z
        mono1 = Monomial([x, y, z], [1, 2, 3])
        mono2 = Monomial([x, y, z], [1, 3, 3])
        mono3 = Monomial([x, z, y, z], [1, 0, 2, 3])

        @test hash(mono1) != hash(mono2)
        @test hash(mono1) == hash(mono3)

        mono4 = Monomial([x], [0])
        @test hash(mono4) == hash(Monomial([], []))
    end
end
