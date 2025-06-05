using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: Polynomial

@testset "Polynomial" begin
	@ncpolyvar x[1:5]
    coeffs = [1.0, 0.0, 2.0, -3.0]
    monos = [Monomial(x, [2, 0, 1, 0, 0]), 
			 Monomial(x, [0, 0, 1, 0, 0]),
			 Monomial(x, [0, 1, 0, 0, 0]),
			 Monomial(x, [0, 0, 0, 1, 0])]

	p = Polynomial(coeffs, monos)

	@test p.coeffs == [2.0, -3.0, 1.0] 
	@test p.monos == monos[[3, 4, 1]] # only non-zero coefficients
end