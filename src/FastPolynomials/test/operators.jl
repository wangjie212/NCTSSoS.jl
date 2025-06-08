using Test, NCTSSoS.FastPolynomials

@testset "^" begin
    @ncpolyvar x y z

    @test x^2 == Monomial([x], [2])
    @test x^0 == Monomial([], [])

    @test_throws AssertionError Base.:(^)(x, -1)
end

@testset "*" begin
    @ncpolyvar x y z

    p1 = 1.0 * x^2
    @test p1 isa Polynomial{Float64}

    p2 = x^2 + y^3
    @test p2 isa Polynomial{Float64}
    @test p2 ≈ Polynomial([1.0, 1.0], [Monomial([x], [2]), Monomial([y], [3])])

    m1 = x * y
    @test m1.vars == [x, y]
    @test m1.z == [1,1]

    m2 = x * y * z
    @test m2.vars == [x, y, z]
    @test m2.z == [1, 1, 1]
end

