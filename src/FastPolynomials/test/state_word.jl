using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials:
    StateWord,
    NCStateWord,
    _projective,
    _unipotent,
    get_state_basis,
    expval,
    neat_dot,
    StatePolynomial

@testset "StateWord and NCStateWord" begin
    @ncpolyvar x[1:2] y[1:2]
    @testset "StateWord" begin
        sw = StateWord([x[1] * x[2], x[2]^2])
        sw2 = ς(x[1] * x[2]) * ς(x[2]^2)
        @test sw == sw2
        @test string(sw) == "<x[1]¹x[2]¹> * <x[2]²>"
        @test sort(variables(sw)) == sort(x)

        sw_sorted = StateWord([x[2]^2, x[1] * x[2]])
        @test sw_sorted == sw

        sw_less = StateWord([x[1] * x[2], x[1]^2])
        @test sw_less < sw

        @test unique([sw, sw_sorted]) == [sw]
        @test sw * sw_less == StateWord([x[2]^2, x[1] * x[2], x[1] * x[2], x[1]^2])

        sw_partial = StateWord([x[1] * x[2]])
        @test sw != sw_partial

        @test degree(sw) == 4

        @test one(sw) == ς(monomial([], []))

        sw_rep1s = StateWord(monomial.(fill(one(x[1]), 3)))
        @test sw_rep1s == StateWord([monomial([], [])])
        @test one(StateWord) == ς(monomial([], []))
        @test (4.0 * sw) isa StatePolynomial{Float64}
    end

    @testset "NCStateWord" begin
        ncsw = NCStateWord(StateWord([x[1] * x[2], x[2]^2]), one(x[1]))
        sw2 = ς(x[1] * x[2]) * ς(x[2]^2)
        @test ncsw.sw == sw2
        @test string(ncsw) == "<x[1]¹x[2]¹> * <x[2]²> * 1"
        @test sort(variables(ncsw)) == sort(x)

        ncsw_sorted = NCStateWord(StateWord([x[2]^2, x[1] * x[2]]), one(x[1]))
        @test ncsw_sorted == ncsw

        ncsw_less = NCStateWord(StateWord([x[1] * x[2], x[1]^2]), one(x[1]))
        @test ncsw_less < ncsw

        @test unique([ncsw, ncsw_sorted]) == [ncsw]
        @test ncsw * ncsw_less ==
            NCStateWord(StateWord([x[2]^2, x[1] * x[2], x[1] * x[2], x[1]^2]), one(x[1]))

        ncsw_partial = NCStateWord(StateWord([x[1] * x[2]]), one(x[1]))
        @test ncsw != ncsw_partial

        @test degree(ncsw) == 4

        @test one(ncsw) == NCStateWord(one(StateWord), one(x[1]))

        ncsw_rep1s = NCStateWord(StateWord(fill(one(Monomial), 3)), one(x[1]))
        @test ncsw_rep1s == NCStateWord(StateWord([one(x[1])]), one(x[1]))

        ncsw1 = NCStateWord(StateWord([x[1] * x[2], x[2]^2]), x[1] * x[2])
        ncsw2 = NCStateWord(StateWord([x[1]^2, x[1]^3]), x[1]^2)
        # NOTE: currently taking <xy>' !=  <xy>
        @test ncsw1' == NCStateWord(StateWord([x[2] * x[1], x[2]^2]), x[2] * x[1])
        @test ncsw1' * ncsw2 == NCStateWord(
            StateWord([x[2] * x[1], x[2]^2, x[1]^2, x[1]^3]), x[2] * x[1] * x[1]^2
        )
        @test neat_dot(ncsw1, ncsw2) ==
            NCStateWord([x[2] * x[1], x[2]^2, x[1]^2, x[1]^3], x[2] * x[1]^3)

        @test expval(ncsw1) == StateWord([x[1] * x[2], x[2]^2, x[1] * x[2]])

            basis = get_state_basis(x, 1, identity)
            total_basis = sort(unique([neat_dot(a, b) for a in basis for b in basis]))
            c_words = [
                [one(x[1])],
                [x[2]],
                [x[2], x[2]],
                [x[2], x[1]],
                [x[1]],
                [x[1], x[1]],
                [one(x[1])],
                [x[2]],
                [x[1]],
                [one(x[1])],
                [x[2]],
                [x[1]],
                [one(x[1])],
                [one(x[1])],
                [one(x[1])],
                [one(x[1])],
            ]
            nc_words =
                monomial.(
                    [
                        fill(one(x[1]), 6)
                        fill(x[2], 3)
                        fill(x[1], 3)
                        [x[2] * x[1], x[2]^2, x[1] * x[2], x[1]^2]
                    ],
                )
            @test sort(total_basis) ==
                sort(map(x -> NCStateWord(x[1], x[2]), zip(c_words, nc_words)))
    end

    @testset "_unipotent" begin
        @ncpolyvar x[1:2] y[1:2]

        @test _unipotent(ς(x[1]^2 * y[1]^2)) == StateWord(Monomial[])
        @test _unipotent(ς(x[1]^2 * y[1] * y[2])) == ς(y[1] * y[2])
        @test _unipotent(ς(x[1]^2) * ς(y[1] * y[2])) == ς(y[1] * y[2])
        @test _unipotent(ς(x[1] * x[2]^2 * x[1]) * ς(y[1] * y[2]^2)) == ς(y[1])
    end

    @testset "_projective" begin
        @ncpolyvar x[1:2] y[1:2]

        @test _projective(ς(x[1]^2 * y[1]^2)) == ς(x[1] * y[1])
        @test _projective(ς(x[1] * x[2]^2 * x[1]) * ς(y[1] * y[2]^2)) ==
            ς(x[1] * x[2] * x[1]) * ς(y[1] * y[2])
    end
end
