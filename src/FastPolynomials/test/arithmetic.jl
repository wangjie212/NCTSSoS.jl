using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: Variable

@testset "Arithmetic" begin
    @testset "Monomial multiplication" begin
        @ncpolyvar x y z

        mono1 = Monomial([x, y], [1, 2])
        mono2 = Monomial([x, z], [3, 4])

        @test mono1 * mono2 == Monomial([x, y, x, z], [1, 2, 3, 4])

        mono3 = Monomial(Variable[], Int64[])
        @test mono1 * mono3 == mono1
        @test mono3 * mono1 == mono1

        mono4 = Monomial([x, y], [1, 0])
        @test mono4 * mono1 == Monomial([x, y], [2, 2])
    end
end
