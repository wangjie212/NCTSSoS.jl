using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials:
    StateWord,
    NCStateWord,
    _projective,
    _unipotent,
    get_state_basis,
    expval,
    neat_dot,
    StatePolynomial,
    Arbitrary, 
    MaxEntangled

const SWord = StateWord{Arbitrary}
const NCSWord = NCStateWord{Arbitrary}
const TWord = StateWord{MaxEntangled}
const NCTWord = NCStateWord{MaxEntangled}

@testset "StateWord and NCStateWord" begin
    @ncpolyvar x[1:2] y[1:2]
    @testset "StateWord" begin
        sw = SWord([x[1] * x[2], x[2]^2])
        sw2 = ς(x[1] * x[2]) * ς(x[2]^2)
        @test sw == sw2

        @test string(sw) == "<x₁¹x₂¹> * <x₂²>"
        @test sort(variables(sw)) == sort(x)

        sw_sorted = SWord([x[2]^2, x[1] * x[2]])
        @test sw_sorted == sw

        sw_less = SWord([x[1] * x[2], x[1]^2])
        @test sw_less < sw

        @test unique([sw, sw_sorted]) == [sw]
        @test sw * sw_less == SWord([x[2]^2, x[1] * x[2], x[1] * x[2], x[1]^2])

        sw_partial = SWord([x[1] * x[2]])
        @test sw != sw_partial

        @test degree(sw) == 4

        @test one(sw) == ς(monomial([], []))

        sw_rep1s = SWord(monomial.(fill(one(x[1]), 3)))
        @test sw_rep1s == SWord([monomial([], [])])
        @test one(SWord) == ς(monomial([], []))
        @test (4.0 * sw) isa StatePolynomial{Float64, Arbitrary}
    end

    @testset "NCStateWord" begin
        ncsw = NCStateWord(SWord([x[1] * x[2], x[2]^2]), one(x[1]))
        sw2 = ς(x[1] * x[2]) * ς(x[2]^2)
        @test ncsw.sw == sw2
        @test string(ncsw) == "<x₁¹x₂¹> * <x₂²> * 1"
        @test sort(variables(ncsw)) == sort(x)

        ncsw_sorted = NCStateWord(SWord([x[2]^2, x[1] * x[2]]), one(x[1]))
        @test ncsw_sorted == ncsw

        ncsw_less = NCStateWord(SWord([x[1] * x[2], x[1]^2]), one(x[1]))
        @test ncsw_less < ncsw

        @test unique([ncsw, ncsw_sorted]) == [ncsw]
        @test ncsw * ncsw_less ==
            NCStateWord(SWord([x[2]^2, x[1] * x[2], x[1] * x[2], x[1]^2]), one(x[1]))

        ncsw_partial = NCStateWord(SWord([x[1] * x[2]]), one(x[1]))
        @test ncsw != ncsw_partial

        @test degree(ncsw) == 4

        @test one(ncsw) == NCStateWord(one(SWord), one(x[1]))

        ncsw_rep1s = NCStateWord(SWord(fill(one(Monomial), 3)), one(x[1]))
        @test ncsw_rep1s == NCStateWord(SWord([one(x[1])]), one(x[1]))

        ncsw1 = NCStateWord(SWord([x[1] * x[2], x[2]^2]), x[1] * x[2])
        ncsw2 = NCStateWord(SWord([x[1]^2, x[1]^3]), x[1]^2)
        # NOTE: currently taking <xy>' !=  <xy>
        @test ncsw1' == NCStateWord(SWord([x[2] * x[1], x[2]^2]), x[2] * x[1])
        @test ncsw1' * ncsw2 == NCStateWord(
            SWord([x[2] * x[1], x[2]^2, x[1]^2, x[1]^3]), x[2] * x[1] * x[1]^2
        )

        @test neat_dot(ncsw1, ncsw2) ==
            NCStateWord(Arbitrary,[x[2] * x[1], x[2]^2, x[1]^2, x[1]^3], x[2] * x[1]^3)

        @test expval(ncsw1) == SWord([x[1] * x[2], x[2]^2, x[1] * x[2]])

        sa = SimplifyAlgorithm(; comm_gps=[x], is_projective=false, is_unipotent=false)
        basis = get_state_basis(Arbitrary, x, 1, sa)
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
            sort(map(x -> NCStateWord(Arbitrary, x[1], x[2]), zip(c_words, nc_words)))
    end

    @testset "_unipotent" begin
        @ncpolyvar x[1:2] y[1:2]

        @test _unipotent(ς(x[1]^2 * y[1]^2)) == SWord(Monomial[])
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
