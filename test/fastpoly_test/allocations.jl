using Test,NCTSSoS.FastPolynomials

using NCTSSoS.FastPolynomials: _comm!
@testset "_comm! allocations" begin
    @ncpolyvar x[1:3] y[1:3]

    comm_gps = Dict(zip([x; y], [fill(1, 3); fill(2, 3)]))
    mono2 = monomial([], [])
	@test (@allocated _comm!(mono2, comm_gps)) === 0
	mono3  = monomial([x[1],y[1],x[2]],[1,2,3])
	@test (@allocated _comm!(mono3, comm_gps)) === 0
end

using NCTSSoS.FastPolynomials: simplify!, SimplifyAlgorithm

@testset "simplify! allocations" begin
    @ncpolyvar x[1:3] y[1:3]

	sa1 = SimplifyAlgorithm(comm_gps=[x,y], is_unipotent=false, is_projective=false)
	sa2 = SimplifyAlgorithm(comm_gps=[x,y], is_unipotent=true, is_projective=false)
	sa3 = SimplifyAlgorithm(comm_gps=[x,y], is_unipotent=false, is_projective=true)

    mono = monomial([x[1], y[1], x[1], y[2], x[2]], [1, 1, 1, 2, 3])

    @test (@allocated simplify!(mono, sa1)) === 0
    @test (@allocated simplify!(mono, sa2)) === 0
    @test (@allocated simplify!(mono, sa3)) === 0
end