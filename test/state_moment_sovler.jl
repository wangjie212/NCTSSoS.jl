using Test, NCTSSoS
using JuMP, DynamicPolynomials
using NCTSSoS: get_state_basis, neat_dot

@test "Constrain Moment matrix" begin
    @ncpolyvar x[1:2]

    basis = get_state_basis(x,1)
    basis = map(x->NCStateWord(StateWord(x[1]),x[2]),basis)

    sp = StatePolynomial([one(Float64)], [StateWord([x[1],x[2]])])
    ncsp = StatePolynomialOp([sp], monomial.([x[1]]))
    poly = one(ncsp)

    total_basis = sort(unique([neat_dot(a,b) for a in basis for b in basis]))
    for b in total_basis
        @show b
    end

    model = GenericModel{Float64}()
    @variable(model, )

end