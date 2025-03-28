using Test, NCTSSoS
using NCTSSoS: StateWord, StatePolynomial, StatePolynomialOp
using DynamicPolynomials

@testset "NCStatePolynomial Components" begin
    @ncpolyvar x[1:2]

    @testset "StateWord" begin
        sw = StateWord([x[1] * x[2], x[2]^2])
        @test string(sw) == "<x[2]^2> * <x[1]*x[2]>"
        @test sort(variables(sw)) == sort(x)
        @test_throws AssertionError StateWord([one(x[1]), x[2]^2, x[1]])

        sw_sorted = StateWord([x[2]^2, x[1] * x[2]])
        @test sw_sorted == sw

        sw_less = StateWord([x[1]*x[2],x[1]^2])
        @test sw_less > sw

        @test unique([sw,sw_sorted]) == [sw]
        @test sw * sw_less == StateWord([x[2]^2, x[1] * x[2], x[1] * x[2], x[1]^2])
    end

    @testset "StatePolynomial" begin
        sws = StateWord.([[x[1]*x[2],x[2]^2],[x[1]*x[2]],[x[2]^3]])
        sp = StatePolynomial([1.0,2.0,5.0],sws)
        @test string(sp) == "1.0 * <x[2]^2> * <x[1]*x[2]> + 2.0 * <x[1]*x[2]> + 5.0 * <x[2]^3>"
        sws_rep = StateWord.([[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^2, x[1] * x[2]], [x[2]^3]])
        sp_rep = StatePolynomial([0.5, 2.0, 0.5, 5.0], sws_rep)
        @test sp == sp_rep

        sws_diff= StateWord.([[x[1]*x[2],x[2]^2],[x[1]*x[2]],[x[2]^3, x[1]]])
        sp_diff = StatePolynomial([1.0,2.0,5.0],sws_diff)
        @test sp_diff != sp
    end
end

@testset "State Polynomial" begin

    @testset "free state polynomial" begin
        @ncpolyvar x[1:2]
        coeffs = [[1.0, 0.5], [2.0], [3.0]]
        formal_words = [StateWord.(wd) for wd in [[[x[1]*x[2],x[2]^2],[x[2]*x[1],x[1]^2]],[[x[1]*x[2]]],[[x[2]^3]]]]
        sps = [StatePolynomial(coef, fw) for (coef, fw) in zip(coeffs, formal_words)]
        words = [one(x[1]), one(x[2]), one(x[1])]
        spop = StatePolynomialOp(sps, words)
        @test string(spop) == "0.5 * <x[2]*x[1]> * <x[1]^2> + 1.0 * <x[2]^2> * <x[1]*x[2]> + 2.0 * <x[1]*x[2]> + 3.0 * <x[2]^3>"
    end

    @testset "nc state polynomial" begin
        @ncpolyvar x[1:2]
        coeffs = [[1.0], [2.0], [3.0]]
        formal_words = [StateWord.(wd) for wd in  [[[x[1]*x[2],x[2]^2]],[[x[1]*x[2]]],[[x[2]^3]]]]
        sps = [StatePolynomial(coef, fw) for (coef, fw) in zip(coeffs, formal_words)]
        words = [x[1], x[2]*x[2], one(x[1])]
        spop = StatePolynomialOp(sps, words)
        @test string(spop) == "3.0 * <x[2]^3> + 1.0 * <x[2]^2> * <x[1]*x[2]> ⋅ x[1] + 2.0 * <x[1]*x[2]> ⋅ x[2]^2"
    end

    @testset "State Polynomial Arithmetic" begin
        @ncpolyvar x[1:2] 

        coeffs = [[1.0, 2.0], [2.0, 3.0]]
        formal_words = [StateWord.(wd) for wd in [[[x[1] * x[2]], [x[2] * x[1]]], [[x[2]^2], [x[1]^2]]]]
        sps = [StatePolynomial(coef, fw) for (coef, fw) in zip(coeffs, formal_words)]
        sp_prod = StatePolynomial([2.0, 4.0, 3.0, 6.0], StateWord.([[x[1]*x[2],x[2]^2],[x[2]*x[1],x[2]^2],[x[1]*x[2],x[1]^2],[x[2]*x[1],x[1]^2]]))
        @test sps[1] * sps[2] == sp_prod
        @test sps[1] + sps[2] == StatePolynomial([1.0, 2.0, 2.0, 3.0], StateWord.([[x[1] * x[2]], [x[2] * x[1]], [x[2]^2], [x[1]^2]]))
    end

    @testset "variables of State Polynomial" begin
        @ncpolyvar x[1:10]

        sws = map(x->StateWord.(x),[[[x[4]*x[5],x[6]*x[7]]],[[x[8],]],[[x[9],x[10]]]])
        coeffs = [[1.0], [2.0], [3.0]]
        sps = [StatePolynomial(coef, fw) for (coef, fw) in zip(coeffs, sws)]
        words = monomial.([x[1], x[2], x[3]])
        sp = StatePolynomialOp(sps, words)
        @test variables(sp) == sort(x)
    end

    @testset "State Polynomial Op comparison" begin
        @ncpolyvar a b c
        words_order1 = [a * b, b * c, c^2]
        words_order2 = [b * c, c^2, a * b]
        state_word_order1 = StateWord.([[a^2], [b^2], [c^2]])
        state_word_order2 = StateWord.([[b^2], [c^2], [a^2]])
        state_poly_order1 = StatePolynomial.([[1.0], [2.0], [3.0]], [[sw] for sw in state_word_order1])
        state_poly_order2 = StatePolynomial.([[2.0], [3.0], [1.0]], [[sw] for sw in state_word_order2])
        spop_order1 = StatePolynomialOp(state_poly_order1, words_order1)
        spop_order2 = StatePolynomialOp(state_poly_order2, words_order2)
        @test spop_order1 == spop_order2
    end
end