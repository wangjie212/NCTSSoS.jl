using Test, NCTSSoS.FastPolynomials
using BenchmarkTools

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

		@ncpolyvar x[1:100000]

		var_coll = sort([x[i] for i in 1:99999])

		@benchmark searchsortedfirst(var_coll, x[100000])	

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

	end

end