using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: star, symmetric_canonicalize, cyclic_canonicalize, support, neat_dot, get_basis,_comm, _projective, _unipotent, monomials

using NCTSSoS.FastPolynomials: Monomial, Polynomial



@testset "Utilities" begin
    @ncpolyvar x y z

    @testset "Star Operation" begin
        mono1 = Monomial([x, y, z], [2, 0, 1])

        # NOTE: I am assuming all variables are Hermitian
        mono1_star = star(mono1)

        @test mono1_star.vars == [z, x]
        @test mono1_star.z == [1, 2]

        mono2 = Monomial([x, y, z], [0, 0, 0])
        mono2_star = star(mono2)

        @test isempty(mono2_star.vars)
        @test isempty(mono2_star.z) 

        mono3 = Monomial([x, y, z], [1, 1, 1]) 
        mono3_star = star(mono3)
        @test mono3_star.vars == [z, y, x] 
        @test mono3_star.z == [1, 1, 1]
    end

    @testset "Symmetric Canonical Form" begin
        mono1 = Monomial([z, y, x], [1, 1, 2]) 

        mono1_sym = symmetric_canonicalize(mono1)
        @test mono1_sym.vars == [x, y, z]
        @test mono1_sym.z == [2, 1, 1]

        mono2 = Monomial([x, y, z], [2, 1, 1])
        mono2_sym = symmetric_canonicalize(mono2)
        @test mono2_sym.vars == [x, y, z]
        @test mono2_sym.z == [2, 1, 1]

        poly1 = Polynomial([0.1, 0.2, 0.3], [Monomial([x, y, z], [2, 1, 1]), Monomial([z, y, x], [1, 1, 2]), Monomial([], [])])
        poly1_sym = symmetric_canonicalize(poly1)

        @test poly1_sym.coeffs ≈ [0.3, 0.3]
        @test poly1_sym.monos == [Monomial([], []), Monomial([x, y, z], [2, 1, 1])]

        n = 3
        @ncpolyvar a[1:n]
        poly3 = Polynomial([1, -1, -1, 3, -2, 2, -1, -1, 6, 9, 9, -54, 142], [
            Monomial([a[1]], [2]),
            Monomial([a[1], a[2]], [1, 1]),
            Monomial([a[2], a[1]], [1, 1]),
            Monomial([a[2]], [2]),
            Monomial([a[1], a[2], a[1]], [1, 1, 1]),
            Monomial([a[1], a[2], a[1]], [1, 2, 1]),
            Monomial([a[2], a[3]], [1, 1]),
            Monomial([a[3], a[2]], [1, 1]),
            Monomial([a[3]], [2]),
            Monomial([a[2], a[3]], [2, 1]),
            Monomial([a[3], a[2]], [1, 2]),
            Monomial([a[3], a[2], a[3]], [1, 1, 1]),
            Monomial([a[3], a[2], a[3]], [1, 2, 1])
        ])

        supp = sort(map([
            [1, 2, 2, 1],
            [3, 2, 2, 3],
            [1, 2, 1],
            [2, 2, 3],
            [3, 2, 3],
            [1, 1],
            [1, 2],
            [2, 2],
            [2, 3],
            [3, 3],
        ]) do v
            Monomial(a[v], ones(length(v)))
        end)

        coe = [-2, 1, -2, 3, 6, -2, 18, -54, 2, 142]

        poly3_sym = symmetric_canonicalize(poly3)

        @test poly3_sym.coeffs == coe
        @test poly3_sym.monos == supp
    end

    @testset "Cyclic Canonical Form" begin
        @test cyclic_canonicalize(x^2*y^2) == x*y^2*x
        @test cyclic_canonicalize(z*y^2) == y*z*y
        @test cyclic_canonicalize(z*z) == z^2
        @test isone(cyclic_canonicalize(one(x)))

        n = 3
        @ncpolyvar a[1:n]
        f = Polynomial([
            1, -1,-1,3,-2,2,-1,-1,6,9,9,-54,142,5,5,5]
            ,[Monomial([a[1]],[2]),
            Monomial([a[1],a[2]],[1,1]),
            Monomial([a[2],a[1]],[1,1]),
            Monomial([a[2]],[2]),
            Monomial([a[1],a[2],a[1]],[1,1,1]),
            Monomial([a[1],a[2],a[1]],[1,2,1]),
            Monomial([a[2],a[3]],[1,1]),
            Monomial([a[3],a[2]],[1,1]),
            Monomial([a[3]],[2]),
            Monomial([a[2],a[3]],[2,1]),
            Monomial([a[3],a[2]],[1,2]),
            Monomial([a[3],a[2],a[3]],[1,1,1]),
            Monomial([a[3],a[2],a[3]],[1,2,1]),
            Monomial([a[1],a[2]],[2,2]),
            Monomial([a[1],a[2],a[3]],[1,1,1]),
            Monomial([a[3],a[2],a[1]],[1,1,1])]
        )

        f2 = Polynomial(
            [7,142,-2,18,-54,1,-2,3,-2,6,10],
            [Monomial([a[1],a[2],a[1]],[1,2,1]),
            Monomial([a[3],a[2],a[3]],[1,2,1]),
            Monomial([a[1],a[2],a[1]],[1,1,1]),
            Monomial([a[2],a[3]],[2,1]),
            Monomial([a[3],a[2],a[3]],[1,1,1]),
            Monomial([a[1]],[2]),
            Monomial([a[1],a[2]],[1,1]),
            Monomial([a[2]],[2]),
            Monomial([a[2],a[3]],[1,1]),
            Monomial([a[3]],[2]),
            Monomial([a[1],a[2],a[3]],[1,1,1])
            ]
        )
        @test cyclic_canonicalize(f) == cyclic_canonicalize(f2)
    end

    @testset "Get basis" begin
        @ncpolyvar x y z

        monomials_deg2 = monomials([x, y, z], 2)
        @test sort(monomials_deg2) == sort([
            Monomial([x], [2]),
            Monomial([y], [2]),
            Monomial([z], [2]),
            Monomial([x, y], [1, 1]),
            Monomial([x, z], [1, 1]),
            Monomial([y, z], [1, 1]),
            Monomial([z, y], [1, 1]),
            Monomial([z, x], [1, 1]),
            Monomial([y, x], [1, 1]),
        ])

        nc_basis_deg2 = get_basis([x, y, z], 2)

        @test sort(nc_basis_deg2) == sort([one(x), Monomial([x], [1]), Monomial([y], [1]), Monomial([z], [1]), x^2, y^2, z^2, x * y, x * z, y * z, z * x, z * y, y * x])
    end

    @testset "support" begin
        @ncpolyvar x y z
        poly = Polynomial(
            [0.1, 0.2, -0.2, 0.3],
            [
                Monomial([x, y], [2, 1]),
                Monomial([], []),
                Monomial([x, y, x], [1, 1, 1]),
                Monomial([z], [1])
            ])

        @test sort(support(poly, identity)) == sort([one(x), Monomial([z], [1]), x^2 * y, x * y * x])
        @test sort(support(poly, cyclic_canonicalize)) == sort([x * y * x, Monomial([z], [1]), one(x)])
    end

    @testset "neat_dot" begin
        @ncpolyvar x y z
        mono1 = Monomial([x, y], [1, 0])

        mono2 = Monomial([x, y], [1, 1])

        @test neat_dot(mono1, mono2) == Monomial([x, y], [2, 1])
        @test neat_dot(mono2, mono2) == Monomial([y, x, y], [1, 2, 1])
    end

    @testset "_comm" begin
        @ncpolyvar a[1:3]
        @ncpolyvar b[1:3]

        mono = a[1]^2*b[2]^2*a[2]*b[1]^3
        @test prod(_comm(mono, [Set(a), Set(b)])) == a[1]^2 * a[2] * b[2]^2 * b[1]^3
        mono = a[1]^3*a[3]
        @test prod(_comm(mono, [Set(a), Set(b)])) == a[1]^3 * a[3]
    end

    @testset "_projective" begin
        mono = y * x^3 * y*z^3
        @test _projective(mono) == y*x*y*z
    end

    @testset "_unipotent" begin
        mono = z*x*y*z^2*y*x*z
        @test _unipotent(mono) == one(mono)
    end


end
