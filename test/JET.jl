using Test, JET, NCTSSoS, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: get_state_basis, MaxEntangled, monomials, get_basis

@testset "Get State Basis Type Stability" begin
    @ncpolyvar x y

	@report_opt monomials([x,y],Val(2))

    err = @report_opt get_basis([x, y], 2)

	err

	@code_warntype get_basis([x, y], 2)

    sa = SimplifyAlgorithm(; comm_gps=[[x, y]], is_projective=false, is_unipotent=false)

    @code_warntype get_state_basis(MaxEntangled, [x, y], 1, sa)
    @report_opt get_state_basis(MaxEntangled, [x, y], 1, sa)
    @report_call get_state_basis(MaxEntangled, [x, y], 1, sa)

end