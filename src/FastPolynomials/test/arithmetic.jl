using Test, NCTSSoS.FastPolynomials

@testset "Arithmetic" begin
	@ncpolyvar x y z

    mono1 = Monomial([x, y], [1, 2])
    mono2 = Monomial([x, z], [3, 4])

    @test mono1 * mono2 == Monomial([x, y, x, z], [1, 2, 3, 4])

	# Test multiplication with a variable
	@test x * mono1 == Monomial([x, y], [1, 2])
	@test mono1 * x == Monomial([x, y], [1, 2])

	# Test multiplication with a constant
	@test 2 * mono1 == Monomial([x, y], [2, 4])
	@test mono1 * 2 == Monomial([x, y], [2, 4])

end