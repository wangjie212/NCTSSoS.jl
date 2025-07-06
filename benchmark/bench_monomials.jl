module BenchMonomials
using BenchmarkTools
using NCTSSoS.FastPolynomials 
using NCTSSoS.FastPolynomials: monomial, monomials

const SUITE = BenchmarkGroup()

@ncpolyvar x[1:10]
var_vec = [x[1], x[2], x[2], x[1], x[3]]
z_vec = [10, 20, 2, 0, 3]

SUITE["Monomial Creation"] = @benchmarkable monomial(var_vec, z_vec)

SUITE["Basis Creation"] = @benchmarkable monomials(x, Val(3))
end

BenchMonomials.SUITE
