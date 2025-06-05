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

    @testset "Polynomial Addition" begin
        @ncpolyvar x y z
        p1 = Polynomial([1.0, 2.0], [Monomial([x], [1]), Monomial([y], [2])])
        p_sum1 = p1 + Monomial([z],[1])
        @test p_sum1 == Polynomial([1.0, 2.0, 1.0], [Monomial([x], [1]), Monomial([y], [2]), Monomial([z], [1])])
        p2 = Monomial([x],[1]) + Monomial([y],[2])
        @test p2.coeffs == [1.0, 1.0]
        @test p2.monos == [Monomial([x], [1]), Monomial([y], [2])]
    end

    @testset "Scaling Polynomial" begin
        @ncpolyvar x y z
        p1 = Polynomial([1.0, 2.0], [Monomial([x], [1]), Monomial([y], [2])])
        p2 = Float32(2.0) * p1
        @test p2.coeffs == [2.0, 4.0] 
        @test p2.monos == [Monomial([x], [1]), Monomial([y], [2])]
    end
end
