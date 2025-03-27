using Test, NCTSSoS
using NCTSSoS: StateWord
using DynamicPolynomials

@testset "StateWord" begin
    @ncpolyvar x[1:2]
    sw = StateWord([x[1]*x[2], x[2]^2])
    @test string(sw) == "⟨x[1]*x[2]⟩ ⟨x[2]^2⟩"
end

@testset "StateTerm" begin

end

@teestset "StatePolynomial" begin

end

@testset "State Polynomial" begin
    @ncpolyvar x[1:2]
    @testset "free state polynomial" begin
        coeffs = [[1.0, 0.5], [2.0], [3.0]]
        formal_words = [[[x[1]*x[2],x[2]^2],[x[2]*x[1],x[1]^2]],[[x[1]*x[2]]],[[one(x[1])]]]
        words = [one(x[1]), one(x[2]), one(x[1])]
        sp = StatePolynomial(coeffs, formal_words, words)
        @test string(sp) == "(1.0 ⋅ ⟨x[1]*x[2]⟩ ⋅ ⟨x[2]^2⟩ + 0.5 ⋅ ⟨x[2]*x[1]⟩ ⋅ ⟨x[1]^2⟩) ⋅ 1 + (2.0 ⋅ ⟨x[1]*x[2]⟩) ⋅ 1 + (3.0 ⋅ ⟨1⟩) ⋅ 1"
    end
    @testset "nc state polynomial" begin
        coeffs = [[1.0], [2.0], [3.0]]
        formal_words = [[[x[1]*x[2],x[2]^2]],[[x[1]*x[2]]],[[one(x[1])]]]
        words = [x[1], x[2]*x[2], one(x[1])]
        sp = StatePolynomial(coeffs, formal_words, words)
        @test string(sp) == "(1.0 ⋅ ⟨x[1]*x[2]⟩ ⋅ ⟨x[2]^2⟩) ⋅ x[1] + (2.0 ⋅ ⟨x[1]*x[2]⟩) ⋅ x[2]^2 + (3.0 ⋅ ⟨1⟩) ⋅ 1"
    end

    @testset "Constructor from vector of polynomials" begin
        @ncpolyvar x[1:2] y[1:2]
        sp = StatePolynomial([[x[1] * y[2] + x[2] * y[1], x[1] * y[2] + x[2] * y[1]], [x[1] * y[1] - x[2] * y[2], x[1] * y[1] - x[2] * y[2]]], [one(x[1]), one(x[1])])
        @test string(sp) == "(1 ⋅ ⟨x[2]*y[1]⟩ ⋅ ⟨x[2]*y[1]⟩ + 1 ⋅ ⟨x[1]*y[2]⟩ ⋅ ⟨x[2]*y[1]⟩ + 1 ⋅ ⟨x[2]*y[1]⟩ ⋅ ⟨x[1]*y[2]⟩ + 1 ⋅ ⟨x[1]*y[2]⟩ ⋅ ⟨x[1]*y[2]⟩) ⋅ 1 + (1 ⋅ ⟨x[2]*y[2]⟩ ⋅ ⟨x[2]*y[2]⟩ + -1 ⋅ ⟨x[1]*y[1]⟩ ⋅ ⟨x[2]*y[2]⟩ + -1 ⋅ ⟨x[2]*y[2]⟩ ⋅ ⟨x[1]*y[1]⟩ + 1 ⋅ ⟨x[1]*y[1]⟩ ⋅ ⟨x[1]*y[1]⟩) ⋅ 1"
    end

    @testset "Utils" begin
        @ncpolyvar x[1:10]

        formal_polynomials = map(x->polynomial.(x),[[x[4]*x[5],x[6]*x[7]],[x[8],],[x[9],x[10]]])
        words = monomial.([x[1], x[2], x[3]])
        sp = StatePolynomial(formal_polynomials, words)
        @test variables(sp) == sort(x)
    end
end