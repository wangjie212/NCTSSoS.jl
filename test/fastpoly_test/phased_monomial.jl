using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: _merge_phase, _phase_to_num

@testset "Phased Monomial" begin
	@testset "Merge Phases" begin
        @test _merge_phase((true, false), (false, true)) == (true, true)
        @test _merge_phase((true, true), (false, true)) == (false, false)
        @test _merge_phase((false, false), (false, false)) == (false, false)
        @test _merge_phase((true, true), (true, true)) == (true, false)
	end
    @testset "Phase to Number" begin
        @test _phase_to_num((true, false)) == -1
        @test _phase_to_num((false, true)) == im
        @test _phase_to_num((true, true)) == -im
        @test _phase_to_num((false, false)) == 1
    end
end