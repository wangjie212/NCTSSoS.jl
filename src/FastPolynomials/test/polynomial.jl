using Test, NCTSSoS.FastPolynomials
# TODO: after API is stabilized, remove function import here
using NCTSSoS.FastPolynomials: Polynomial, variables

@testset "Polynomial" begin
    @testset "creation" begin
        @ncpolyvar x[1:5]
        coeffs = [1.0, 0.0, 2.0, -3.0]
        monos = [
            Monomial(x, [2, 0, 1, 0, 0]),
            Monomial(x, [0, 0, 1, 0, 0]),
            Monomial(x, [0, 1, 0, 0, 0]),
            Monomial(x, [0, 0, 0, 1, 0]),
        ]

        p = Polynomial(coeffs, monos)

        @test p.coeffs == [2.0, -3.0, 1.0]
        @test p.monos == monos[[3, 4, 1]] # only non-zero coefficients

        p_rep_mono = Polynomial([1.0,2.0],[Monomial([x[1],x[2]],[1,2]),Monomial([x[1],x[2],x[2]],[1,1,1])])

        @test p_rep_mono.coeffs == [3.0]
        @test p_rep_mono.monos == [Monomial([x[1], x[2]], [1, 2])]
    end
    @testset "utils" begin
        @ncpolyvar x y z
        p = Polynomial([1.0,2.0,3.0], [Monomial([x,y], [1,2]), Monomial([y,z], [2,3]), Monomial([z,x], [3,4])])

        @test variables(p) == [x,y,z]
    end
end
