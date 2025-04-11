using Test, NCTSSoS
using NCTSSoS: StateWord, StateTerm, StatePolynomial, NCStateWord, NCStateTerm, StatePolynomialOp, expval, neat_dot, get_state_basis
using DynamicPolynomials
using DynamicPolynomials: monomial

@testset "NCStatePolynomial Components" begin
    @ncpolyvar x[1:2] y[1:2]

    @testset "StateWord" begin
        sw = StateWord([x[1] * x[2], x[2]^2])
        sw2 = ς(x[1]*x[2]) * ς(x[2]^2)
        @test sw == sw2
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
        @test (4.0 * sw) isa StateTerm 
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

        basis = get_state_basis(x,1, identity)
        total_basis = sort(unique([neat_dot(a,b) for a in basis for b in basis]))
        c_words = map(x->StateWord(monomial.(x)),[[one(x[1])],[x[2]],[x[2],x[2]],[x[2],x[1]],[x[1]],[x[1],x[1]],[one(x[1])],[x[2]],[x[1]],[one(x[1])],[x[2]],[x[1]],[one(x[1])],[one(x[1])],[one(x[1])],[one(x[1])]])
        nc_words = monomial.([fill(one(x[1]),6);fill(x[2],3);fill(x[1],3);[x[2]*x[1],x[2]^2,x[1]*x[2],x[1]^2]])
        @test total_basis == map(x->NCStateWord(x[1],x[2]),zip(c_words,nc_words))

        @test 2.0 * ncsw1 isa NCStateTerm
    end

    @testset "StatePolynomial" begin
        sws = StateWord.([[x[1]*x[2],x[2]^2],[x[1]*x[2]],[x[2]^3]])
        sts = StateTerm.([1.0, 2.0, 5.0], sws)
        sp = StatePolynomial(sts)
        @test string(sp) == "1.0 * <x[2]^2> * <x[1]*x[2]> + 2.0 * <x[1]*x[2]> + 5.0 * <x[2]^3>"
        sws_rep = StateWord.([[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^2, x[1] * x[2]], [x[2]^3]])
        sts_rep = StateTerm.([0.5, 2.0, 0.5, 5.0], sws_rep)
        sp_rep = StatePolynomial(sts_rep)
        @test sp == sp_rep

        sws_diff = StateWord.([[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^3, x[1]]])
        sts_diff = StateTerm.([1.0,2.0,5.0],sws_diff)
        sp_diff = StatePolynomial(sts_diff)
        @test sp_diff != sp

        @test degree(sp) == 4

        @test one(sp) == StatePolynomial([StateTerm(1.0, one(sws[1]))])

        @test monomials(sp) == StateWord.([[x[2]^2, x[1] * x[2]], [x[1] * x[2]], [x[2]^3]])

        sp_ez = -1.0 *  ς(x[1]*y[1]) * one(x[1]) - 1.0 * ς(x[1]*y[2]) * one(x[1]) - (1.0 * ς(x[2]*y[1]) * one(x[1]) ) + 1.0 * ς(x[2]*y[2]) * one(x[1])
        sp_hard = StatePolynomialOp(map(a -> a[1]*StateWord([a[2]])*a[3], zip([-1.0, -1.0, -1.0, 1.0], [x[1] * y[1], x[1] * y[2], x[2] * y[1], x[2] * y[2]],fill(one(x[1]),4))))
        @test sp_ez == sp_hard
    end

    @testset "StateTerm" begin
        # TODO: add test if API persists
    end

    @testset  "NCStateTerm" begin
        # TODO: add test if API persists 
    end
end


@testset "State Polynomial" begin
    @testset "free state polynomial" begin
        @ncpolyvar x[1:2]
        sws = [prod(ς.(wd)) for wd in [[x[1]*x[2],x[2]^2],[x[2]*x[1],x[1]^2],[x[1]*x[2]],[x[2]^3]]]
        ncsws = [NCStateWord(sw,monomial(one(x[1]))) for sw in sws]
        spop = StatePolynomialOp([NCStateTerm(coef, ncsw) for (coef, ncsw) in zip([1.0, 0.5, 2.0, 3.0], ncsws)])
        string(spop)
        @test string(spop) == "0.5 * <x[2]*x[1]> * <x[1]^2> ⋅ 1 + 1.0 * <x[2]^2> * <x[1]*x[2]> ⋅ 1 + 2.0 * <x[1]*x[2]> ⋅ 1 + 3.0 * <x[2]^3> ⋅ 1"
    end

    @testset "nc state polynomial" begin
        @ncpolyvar x[1:2]
        formal_words = [prod(ς.(wd)) for wd in  [[x[1]*x[2],x[2]^2],[x[1]*x[2]],[x[2]^3]]]
        ncsws = [sw*ncw for (sw, ncw) in zip(formal_words, [x[1], x[2]^2, one(x[1])])]
        spop = StatePolynomialOp([coef * ncsw for (coef, ncsw) in zip([1.0, 2.0, 3.0], ncsws)])
        @test string(spop) == "3.0 * <x[2]^3> ⋅ 1 + 1.0 * <x[2]^2> * <x[1]*x[2]> ⋅ x[1] + 2.0 * <x[1]*x[2]> ⋅ x[2]^2"

        @test degree(spop) == 5
        @test monomials(spop) == map(a -> prod(ς.(a[1])) * monomial(a[2]), [([x[2]^3], one(x[1])), ([x[2]^2, x[1] * x[2]], x[1]), ([x[1] * x[2]], x[2]^2)])
    end

    @testset "State Polynomial Arithmetic" begin
        @ncpolyvar x[1:2] 

        sws = [StateTerm(coef, StateWord([wd])) for (coef, wd) in zip([1.0, 2.0, 2.0, 3.0], [x[1] * x[2], x[2] * x[1], x[2]^2, x[1]^2])]
        sps = [sws[1]+sws[2],sws[3]+sws[4]]

        sp_prod = StatePolynomial([StateTerm(coef,sw) for (coef,sw) in zip([2.0, 4.0, 3.0, 6.0], StateWord.([[x[1]*x[2],x[2]^2],[x[2]*x[1],x[2]^2],[x[1]*x[2],x[1]^2],[x[2]*x[1],x[1]^2]]))])
        @test sps[1] * sps[2] == sp_prod
        sp_sum = StatePolynomial([StateTerm(coef,sw) for (coef,sw) in zip([1.0, 2.0, 2.0, 3.0], StateWord.([[x[1] * x[2]], [x[2] * x[1]], [x[2]^2], [x[1]^2]]))])
        @test sps[1] + sps[2] == sp_sum
    end

    @testset "variables of State Polynomial Op" begin
        @ncpolyvar x[1:10]

        sts = map(x->StateTerm(x[1],StateWord(x[2])),zip([1.0,2.0,3.0],[[x[4]*x[5],x[6]*x[7]],[x[8]],[x[9],x[10]]]))
        sp = StatePolynomial(sts)
        @test sort(variables(sp)) == sort(x[4:10])
        words = monomial.([x[1], x[2], x[3]])
        spop = StatePolynomialOp(NCStateTerm.([1.0,2.0,3.0],NCStateWord.(StateWord.([[x[4]*x[5],x[6]*x[7]],[x[8]],[x[9],x[10]]]),words)))
        @test variables(spop) == sort(x)
    end

    @testset "State Polynomial Op comparison" begin
        @ncpolyvar a b c
        words_order1 = [a * b, b * c, c^2]
        words_order2 = [b * c, c^2, a * b]
        state_word_order1 = StateWord.([[a^2], [b^2], [c^2]])
        state_word_order2 = StateWord.([[b^2], [c^2], [a^2]])
        state_poly_order1 = StatePolynomial([StateTerm(coef,sw) for (coef,sw) in zip([1.0, 2.0, 3.0], state_word_order1)])
        state_poly_order2 = StatePolynomial([StateTerm(coef,sw) for (coef,sw) in zip([2.0, 3.0,1.0],  state_word_order2)])
        @test state_poly_order1 == state_poly_order2
        spop_order1 = StatePolynomialOp([NCStateTerm(coef,NCStateWord(sw,monomial(one(a)))) for (coef,sw) in zip([1.0, 2.0, 3.0], state_word_order1)])
        spop_order2 = StatePolynomialOp([NCStateTerm(coef,NCStateWord(sw,monomial(one(a)))) for (coef,sw) in zip([2.0, 3.0,1.0],  state_word_order2)]) 
        @test spop_order1 == spop_order2
    end
end

@testset "Utils" begin
    @ncpolyvar x y 

    get_state_basis([x, y], 1, identity)
    c_words = [[one(x)], [y], [x], [one(x)], [one(x)]]
    nc_words = [one(x), one(x), one(x), y, x]
    @test sort(get_state_basis([x, y], 1, identity)) == sort(map(x -> StateWord(x[1]) * x[2], zip(c_words, nc_words)))
    c_words = [[one(x)], [y], [x], [y * x], [y^2], [x^2], [y, y], [y, x], [x, x], [one(x)], [y], [x], [one(x)], [y], [x], fill([one(x)], 4)...]
    nc_words = [fill(one(x), 9); fill(y, 3); fill(x, 3); [y * x, y^2, x * y, x^2]]
    @test sort(get_state_basis([x, y], 2, identity)) == sort(map(x -> StateWord(x[1]) * x[2], zip(c_words, nc_words)))

    nc_words = [fill(one(x),7);fill(x,4);fill(x^2,2);[x^3]]
    c_words = [[one(x)], [x], [x^2], [x^3], [x, x], [x, x^2], [x, x, x], [one(x)], [x], [x^2], [x, x], [one(x)], [x], [one(x)]]
    @test sort(get_state_basis([x], 3, identity)) == sort(map(x->StateWord(x[1])*x[2],zip(c_words,nc_words))) 
end