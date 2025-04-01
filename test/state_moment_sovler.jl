using Test, NCTSSoS
using JuMP, DynamicPolynomials
using NCTSSoS: get_state_basis, neat_dot, NCStateWord, StateWord, StatePolynomial, StatePolynomialOp, constrain_moment_matrix!, expval, substitute_variables

@testset "Constrain Moment matrix" begin
    @ncpolyvar x[1:2]

    basis = get_state_basis(x,1)
    basis = map(x->NCStateWord(StateWord(x[1]),x[2]),basis)

    sp = map(a -> StatePolynomial([a[1]], [a[2]]), zip([1.0, 2.0, 3.0], StateWord.([[x[1], x[2]], [x[1]], [x[2]]])))
    nc_words = monomial.([one(x[1]), x[1], x[2]])
    ncsp = StatePolynomialOp(sp, nc_words)
    poly = one(ncsp)

    total_basis = sort(unique([expval(neat_dot(a,b)) for a in basis for b in basis]))

    model = GenericModel{Float64}()
    @variable(model, y[1:length(total_basis)])
    wordmap = Dict(zip(total_basis,y))

    ncwords = map(a -> NCStateWord(terms(a[1])[1][2], a[2]), zip(sp, nc_words))
    @test sort(terms(ncsp)) == sort(collect(zip([1.0, 2.0, 3.0], ncwords)))

    @test substitute_variables(ncsp, wordmap) == 1.0 * y[4] + 3.0 * y[3] + 2.0 * y[6]

    true_mom_mtx = expval.([neat_dot(a,b) for a in basis, b in basis])
    mom_mtx_cons = constrain_moment_matrix!(model, one(ncsp), basis, wordmap, PSDCone(), identity)
    mom_mtx = constraint_object(mom_mtx_cons)
    @test reshape(mom_mtx.func, 5, 5) == AffExpr[y[1] y[2] y[5] y[2] y[5]; y[2] y[3] y[4] y[3] y[4]; y[5] y[4] y[6] y[4] y[6]; y[2] y[3] y[4] y[8] y[7]; y[5] y[4] y[6] y[9] y[10]]
end