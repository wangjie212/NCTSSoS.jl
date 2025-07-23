using Test,NCTSSoS.FastPolynomials

using NCTSSoS.FastPolynomials: _comm!
@test "_comm! allocations" begin
    @ncpolyvar x[1:3] y[1:3]

    comm_gps = Dict(zip([x; y], [fill(1, 3); fill(2, 3)]))
    mono = monomial([x[1]], [1])
    @test (@allocated _comm!(mono, comm_gps)) == 0
    mono2 = monomial([], [])
    @test (@allocated _comm!(mono, comm_gps)) == 0
	mono3  = monomial([x[1],y[1],x[2]],[1,2,3])
    @test (@allocated _comm!(mono, comm_gps)) == 0
end