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

SUITE["Compare different degree"] = @benchmarkable cmp(a, b) setup = (a = monomial([x[1], x[2]], [2, 1]); b = monomial([x[1], x[2]], [1, 1])) # TODO: each 1.5ns

SUITE["Compare same degree"] = @benchmarkable cmp(a, b) setup = (a = monomial([x[1], x[2],x[1]], [2, 1, 1]); b = monomial([x[1], x[2]], [3, 1])) # TODO: reach 2.38ns
end

BenchMonomials.SUITE
