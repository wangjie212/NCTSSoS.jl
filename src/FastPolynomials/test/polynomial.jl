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
        @test (-p).coeffs == [-2.0, 3.0, -1.0]
        @test (-p).monos == monos[[3, 4, 1]] # negation preserves monomials


        p_rep_mono = Polynomial([1.0,2.0],[Monomial([x[1],x[2]],[1,2]),Monomial([x[1],x[2],x[2]],[1,1,1])])

        @test p_rep_mono.coeffs == [3.0]
        @test p_rep_mono.monos == [Monomial([x[1], x[2]], [1, 2])]

        p_const = Polynomial(Float32, 1)
        @test p_const isa Polynomial{Float32}
        @test p_const.coeffs == [1.0]
        @test p_const.monos == [Monomial([], [])]
    end

    @testset "creation2" begin
        order = 3
        n = 2
        @ncpolyvar x[1:n]
        f = 0.0
        for i in 1:n
            jset = max(1, i - 5):min(n, i + 1)
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
        coeffs = [25.0, 25.0, -10.0, -10.0, 21.0, -10.0, 21.0, -10.0, 12.0, -4.0, 12.0, -4.0, 3.0, -4.0, 3.0, -4.0, 2.0, 2.0, 2.0]
        mono_zs = [[6, 0, 0], [0, 6, 0], [3, 2, 0], [0, 3, 2], [4, 0, 0], [3, 1, 0], [0, 4, 0], [0, 3, 1], [3, 0, 0], [1, 2, 0], [0, 3, 0], [0, 1, 2], [2, 0, 0], [1, 1, 0], [0, 2, 0], [0, 1, 1], [1, 0, 0], [0, 1, 0], [0, 0, 0]]
        mono_vars = [[1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1], [1, 2, 1]]

        f_monos = [
            Monomial(x[var_idcs], var_zs) for (var_idcs, var_zs) in zip(mono_vars, mono_zs)
        ]
        f_dynamic = mapreduce(+,zip(coeffs,f_monos)) do (coef, mono)
            coef * mono
        end
        @test f_dynamic == f 
    end

    @testset "exponentiate" begin
        @ncpolyvar x
        p = (2x + 5x^3 +1)^2
        @test p == (1 + 4x + 4x^2 +10x^3 + 20x^4 + 25x^6)
    end
    @testset "Addition" begin
        @ncpolyvar x y
        p = 0.0
        p = p + x
        @test p == Polynomial(x)
    end
    @testset "Promotion" begin
        @ncpolyvar x
        @test Polynomial(0.0) == Polynomial([0.0],[Monomial([x], [0])])
    end
    @testset "utils" begin
        @ncpolyvar x y z
        p = Polynomial([1.0,2.0,3.0], [Monomial([x,y], [1,2]), Monomial([y,z], [2,3]), Monomial([z,x], [3,4])])
        @test variables(p) == [x,y,z]
    end
end
