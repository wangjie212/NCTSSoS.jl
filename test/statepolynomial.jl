using Test, NCTSSoS
using DynamicPolynomials

@testset "State Polynomial" begin
    @ncpolyvar x[1:2]
    @testset "free state polynomial" begin
        coeffs = [1.0, 2.0, 3.0]
        formal_words = [[x[1]*x[2],x[2]^2],[x[1]*x[2]],[one(x[1])]]
        words = [one(x[1]), one(x[2]), one(x[1])]
        sp = StatePolynomial(coeffs, formal_words, words)
        @test string(sp) == "1.0 ⋅ ⟨x[1]*x[2]⟩ ⋅ ⟨x[2]^2⟩ ⋅ 1 + 2.0 ⋅ ⟨x[1]*x[2]⟩ ⋅ 1 + 3.0 ⋅ ⟨1⟩ ⋅ 1"
    end
    @testset "nc state polynomial" begin
        coeffs = [1.0, 2.0, 3.0]
        formal_words = [[x[1]*x[2],x[2]^2],[x[1]*x[2]],[one(x[1])]]
        words = [x[1], x[2]*x[2], one(x[1])]
        sp = StatePolynomial(coeffs, formal_words, words)
        @test string(sp) == "1.0 ⋅ ⟨x[1]*x[2]⟩ ⋅ ⟨x[2]^2⟩ ⋅ x[1] + 2.0 ⋅ ⟨x[1]*x[2]⟩ ⋅ x[2]^2 + 3.0 ⋅ ⟨1⟩ ⋅ 1"
    end
end