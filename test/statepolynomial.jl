using Test, NCTSSoS
using NCTSSoS: StateWord, StatePolynomial, StatePolynomialOp, NCStateWord, expval, neat_dot
using DynamicPolynomials
using DynamicPolynomials: monomial

@testset "NCStatePolynomial Components" begin
    @ncpolyvar x[1:2]

    @testset "StateWord" begin
        sw = StateWord([x[1] * x[2], x[2]^2])
        @test string(sw) == "<x[2]^2> * <x[1]*x[2]>"
        @test sort(variables(sw)) == sort(x)

        sw_sorted = StateWord([x[2]^2, x[1] * x[2]])
        @test sw_sorted == sw

        sw_less = StateWord([x[1]*x[2],x[1]^2])
        @test sw_less > sw

        @test unique([sw,sw_sorted]) == [sw]
        @test sw * sw_less == StateWord([x[2]^2, x[1] * x[2], x[1] * x[2], x[1]^2])

        sw_partial = StateWord([x[1]*x[2]])
        @test sw != sw_partial

        @test degree(sw) == 4

        @test one(sw) == StateWord(monomial.([one(x[1])]))

        sw_rep1s = StateWord(monomial.(fill(one(x[1]),3)))
        @test sw_rep1s == StateWord([monomial(one(x[1]))])
    end

    @testset "NCStateWord" begin
        sw1 = StateWord([x[1] * x[2], x[2]^2])
        sw2 = StateWord([x[1]^2, x[1]^3])
        ncsw1 = NCStateWord(sw1, x[1] * x[2])
        ncsw2 = NCStateWord(sw2, x[1]^2)
        @test ncsw1' == NCStateWord(sw1, x[2]*x[1])
        @test ncsw1' * ncsw2 == NCStateWord(sw1 * sw2, x[2] * x[1] * x[1]^2)
        @test neat_dot(ncsw1, ncsw2) == NCStateWord(sw1 * sw2, x[2] * x[1] * x[1]^2)

        @test expval(ncsw1) == StateWord([x[1]*x[2],x[2]^2, x[1]*x[2]])


        @ncpolyvar x[1:2]

        basis = get_state_basis(x,1)
        basis = map(x->NCStateWord(StateWord(x[1]),x[2]),basis)
        total_basis = sort(unique([neat_dot(a,b) for a in basis for b in basis]))
        c_words = map(x->StateWord(monomial.(x)),[[one(x[1])],[x[2]],[x[2],x[2]],[x[2],x[1]],[x[1]],[x[1],x[1]],[one(x[1])],[x[2]],[x[1]],[one(x[1])],[x[2]],[x[1]],[one(x[1])],[one(x[1])],[one(x[1])],[one(x[1])]])
        nc_words = monomial.([fill(one(x[1]),6);fill(x[2],3);fill(x[1],3);[x[2]*x[1],x[2]^2,x[1]*x[2],x[1]^2]])
        @test total_basis == map(x->NCStateWord(x[1],x[2]),zip(c_words,nc_words))
        total_basis
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

        @test degree(sp) == 4

        @test one(sp) == StatePolynomial([one(Float64)],[one(sw)])

        sp

        @test monomials(sp) == StateWord.([[x[2]^2, x[1] * x[2]], [x[1] * x[2]], [x[2]^3]])
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

        @test degree(spop) == 5
        @test monomials(spop) == map(a->NCStateWord(StateWord(a[1]),monomial(a[2])),[([x[2]^3],one(x[1])),([x[2]^2,x[1]*x[2]],x[1]),([x[1]*x[2]],x[2]^2)])
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

@testset "Utils" begin
    @ncpolyvar x y 

    NCTSSoS.get_state_basis([x, y], 1)
    c_words = [[one(x)], [y], [x], [one(x)], [one(x)]]
    nc_words = [one(x), one(x), one(x),y,x]
    @test sort(NCTSSoS.get_state_basis([x, y], 1)) == sort(map(x->(monomial.(x[1]),x[2]),zip(c_words,nc_words))) 
    c_words = [[one(x)], [y], [x], [y * x], [y^2], [x * y], [x^2], [y, y], [y, x], [x, x], [one(x)],[y],[x],[one(x)],[y],[x],fill([one(x)],4)...]
    nc_words = [fill(one(x),10);fill(y,3);fill(x,3);[y*x,y^2,x*y,x^2]]
    @test sort(NCTSSoS.get_state_basis([x, y], 2)) == sort(map(x->(monomial.(x[1]),x[2]),zip(c_words,nc_words))) 

    nc_words = [fill(one(x),7);fill(x,4);fill(x^2,2);[x^3]]
    c_words = [[one(x)], [x], [x^2], [x^3], [x, x], [x, x^2], [x, x, x], [one(x)], [x], [x^2], [x, x], [one(x)], [x], [one(x)]]
    @test sort(NCTSSoS.get_state_basis([x], 3)) == sort(map(x->(monomial.(x[1]),x[2]),zip(c_words,nc_words))) 
end