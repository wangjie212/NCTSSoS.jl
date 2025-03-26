using Test, NCTSSoS

@testset "State Polynomial" begin
    @ncpolyvar x[1:2]
    @testset "Example 7.2.1" begin
        coeffs = [1.0, 2.0, 3.0]
        formal_words = [[x[1]*x[2],x[2]*x[3]],[x[1]*x[2],x[2]*x[3]],[x[1]*x[2],x[2]*x[3]]]
    end
end