using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials:
    StateWord,
    NCStateWord,
    degree,
    StatePolynomial,
    get_state_basis,
    neat_dot,
    expval,
    monomials,
    _unipotent,
    _projective

@testset "NCStatePolynomial Components" begin
    @ncpolyvar x[1:2] y[1:2]

    @testset "NCStatePolynomial" begin
        sws = NCStateWord.([[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^3]], Ref(one(x[1])))
        sp = ncstatepoly([1.0, 2.0, 5.0], sws)
        @test string(sp) ==
            "2.0 * <x[1]¹x[2]¹> * 1 + 5.0 * <x[2]³> * 1 + 1.0 * <x[1]¹x[2]¹> * <x[2]²> * 1"
        sws_rep = NCStateWord.(
            [[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^2, x[1] * x[2]], [x[2]^3]],
            Ref(one(x[1])),
        )
        sp_rep = ncstatepoly([0.5, 2.0, 0.5, 5.0], sws_rep)
        @test sp == sp_rep

        sws_diff = NCStateWord.(
            [[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^3, x[1]]], Ref(one(x[1]))
        )
        sp_diff = ncstatepoly([1.0, 2.0, 5.0], sws_diff)
        @test sp_diff != sp

        @test degree(sp) == 4

        @test one(sp) == ncstatepoly([1.0], [one(sws[1])])

        @test sort(monomials(sp)) == sort(
            NCStateWord.([[x[2]^2, x[1] * x[2]], [x[1] * x[2]], [x[2]^3]], Ref(one(x[1])))
        )

        sp_ez =
            -1.0 * ς(x[1] * y[1]) * one(x[1]) - 1.0 * ς(x[1] * y[2]) * one(x[1]) -
            (1.0 * ς(x[2] * y[1]) * one(x[1])) + 1.0 * ς(x[2] * y[2]) * one(x[1])
        sp_hard = mapreduce(
            a -> a[1] * NCStateWord([a[2]], a[3]),
            +,
            zip(
                [-1.0, -1.0, -1.0, 1.0],
                [x[1] * y[1], x[1] * y[2], x[2] * y[1], x[2] * y[2]],
                fill(one(x[1]), 4),
            ),
        )
        @test sp_ez == sp_hard
    end
end

@testset "State Polynomial" begin
    @testset "free state polynomial" begin
        @ncpolyvar x[1:2]
        ncsws = [
            prod(ς.(wd)) for
            wd in [[x[1] * x[2], x[2]^2], [x[2] * x[1], x[1]^2], [x[1] * x[2]], [x[2]^3]]
        ]
        spop = StatePolynomial([1.0, 0.5, 2.0, 3.0], ncsws)
        @test string(spop) ==
            "2.0 * <x[1]¹x[2]¹> + 3.0 * <x[2]³> + 1.0 * <x[1]¹x[2]¹> * <x[2]²> + 0.5 * <x[1]²> * <x[2]¹x[1]¹>"
    end

    @testset "nc state polynomial" begin
        @ncpolyvar x[1:2]
        ncsws = [
            NCStateWord(wd, ncw) for (wd, ncw) in
            zip([[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^3]], [x[1], x[2]^2, one(x[1])])
        ]

        spop = ncstatepoly([1.0, 2.0, 3.0], ncsws)

        @test string(spop) ==
            "3.0 * <x[2]³> * 1 + 2.0 * <x[1]¹x[2]¹> * x[2]² + 1.0 * <x[1]¹x[2]¹> * <x[2]²> * x[1]¹"

        @test degree(spop) == 5
        @test sort(monomials(spop)) == sort(
            map(
                a -> NCStateWord(a[1], monomial(a[2])),
                [
                    ([x[2]^3], one(x[1])),
                    ([x[2]^2, x[1] * x[2]], x[1]),
                    ([x[1] * x[2]], x[2]^2),
                ],
            ),
        )
    end

    @testset "State Polynomial Arithmetic" begin
        @ncpolyvar x[1:2]

        sps = [
            1.0 * ς(x[1] * x[2]) + 2.0 * ς(x[2] * x[1]), 2.0 * ς(x[2]^2) + 3.0 * ς(x[1]^2)
        ]

        sp_prod =
            2.0 * ς(x[1] * x[2]) * ς(x[2]^2) +
            4.0 * ς(x[2] * x[1]) * ς(x[2]^2) +
            3.0 * ς(x[1] * x[2]) * ς(x[1]^2) +
            6.0 * ς(x[2] * x[1]) * ς(x[1]^2)
        @test sps[1] * sps[2] == sp_prod
        sp_sum =
            1.0 * ς(x[1] * x[2]) + 2.0 * ς(x[2] * x[1]) + 2.0 * ς(x[2]^2) + 3.0 * ς(x[1]^2)
        @test sps[1] + sps[2] == sp_sum
    end

    @testset "variables of NC State Polynomial" begin
        @ncpolyvar x[1:10]

        sp = 1.0 * (ς(x[4] * x[5]) * ς(x[6] * x[7])) + 2.0 * ς(x[8]) + 3.0 * ς(x[9] * x[10])
        @test sort(variables(sp)) == sort(x[4:10])

        words = monomial.([x[1], x[2], x[3]])
        spop =
            1.0 * NCStateWord([x[4] * x[5], x[6] * x[7]], x[1]) +
            2.0 * NCStateWord([x[8]], x[2]) +
            3.0 * NCStateWord([x[9], x[10]], x[3])
        @test variables(spop) == sort(x)
    end

    @testset "State Polynomial Op comparison" begin
        @ncpolyvar a b c
        words_order1 = [a * b, b * c, c^2]
        words_order2 = [b * c, c^2, a * b]
        state_word_order1 = NCStateWord.([[a^2], [b^2], [c^2]], Ref(one(a)))
        state_word_order2 = NCStateWord.([[b^2], [c^2], [a^2]], Ref(one(a)))
        state_poly_order1 = ncstatepoly([1.0, 2.0, 3.0], state_word_order1)
        state_poly_order2 = ncstatepoly([2.0, 3.0, 1.0], state_word_order2)
        @test state_poly_order1 == state_poly_order2
    end
end

@testset "get state basis" begin
    @ncpolyvar x y

    sa = SimplifyAlgorithm(; comm_gps=[[x, y]], is_projective=false, is_unipotent=false)
    get_state_basis([x, y], 1, sa)
    c_words = [[one(x)], [y], [x], [one(x)], [one(x)]]
    nc_words = [one(x), one(x), one(x), y, x]
    @test sort(get_state_basis([x, y], 1, sa)) ==
        sort(map(x -> NCStateWord(x[1], x[2]), zip(c_words, nc_words)))

    c_words = [
        [one(x)],
        [y],
        [x],
        [x * y],
        [y^2],
        [x^2],
        [y, y],
        [y, x],
        [x, x],
        [one(x)],
        [y],
        [x],
        [one(x)],
        [y],
        [x],
        fill([one(x)], 4)...,
    ]
    nc_words = [fill(one(x), 9); fill(y, 3); fill(x, 3); [y * x, y^2, x * y, x^2]]
    @test sort(get_state_basis([x, y], 2, sa)) ==
        sort(map(x -> NCStateWord(x[1], x[2]), zip(c_words, nc_words)))

    nc_words = [fill(one(x), 7); fill(x, 4); fill(x^2, 2); [x^3]]
    c_words = [
        [one(x)],
        [x],
        [x^2],
        [x^3],
        [x, x],
        [x, x^2],
        [x, x, x],
        [one(x)],
        [x],
        [x^2],
        [x, x],
        [one(x)],
        [x],
        [one(x)],
    ]
    @test sort(get_state_basis([x], 3, sa)) ==
        sort(map(x -> NCStateWord(x[1], x[2]), zip(c_words, nc_words)))
end
