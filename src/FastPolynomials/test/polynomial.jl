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

	@test p.coeffs == [1.0, 2.0, -3.0] 
	@test p.monos == monos[[1, 3, 4]] # only non-zero coefficients
end