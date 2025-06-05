using Test, NCTSSoS.FastPolynomials

@testset "Comparison" begin
	@testset "cmp variables" begin
		@ncpolyvar x y z

        @test cmp(x, x) == 0
		@test cmp(x, y) == -1
		@test cmp(y, x) == 1 

		@test isless(x, y) == true
		@test sort([z,y,x]) == [x,y,z]

		@ncpolyvar x[1:10]
		@test cmp(x[1],x[10]) == 1


		@test x[1] == x[1]
		@test x[1] != x[2]
	end

	@testset "cmp monomials" begin
		@ncpolyvar x y z

		mono1 = Monomial([x, y], [1, 2]) 
		mono2 = Monomial([x, y], [1, 2])
		mono3 = Monomial([x, y], [1, 1])
       	mono4 = Monomial([x, z, y], [1, 0, 2])
        mono5 = Monomial([x, z, y], [1, 1, 1])

		# TODO: add more tests

		@test cmp(mono1, mono2) == 0	
		@test cmp(mono1, mono3) == 1		
		@test cmp(mono1,mono4) == 0
		@test cmp(mono4, mono5) == -1

		@test isless(mono4, mono5)
	end

	@testset "compare polynomials" begin
		p1 = Polynomial([1.0 ,2.0], [Monomial([x, y], [1, 2]), Monomial([x, y], [1, 1])])
		p2 = Polynomial([1.0 ,2.0], [Monomial([x, y], [1, 2]), Monomial([x, y], [1, 1])])
		p3 = Polynomial([1.00000001 ,2.0], [Monomial([x, y], [1, 2]), Monomial([x, y], [1, 1])])
		@test p1 == p2
		@test !(p1 == p3)

		@test !isapprox(p1,p3;atol=1e-9)
		@test isapprox(p1,p3;atol=1e-7)
	end
end