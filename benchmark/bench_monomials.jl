module BenchMonomials
using BenchmarkTools
using NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: monomials

@ncpolyvar x[1:10]

const SUITE = BenchmarkGroup()

var_vec = [x[1], x[2], x[2], x[1], x[3]]
z_vec = [10, 20, 2, 0, 3]

SUITE["Monomial Creation"] = @benchmarkable Monomial(var_vec, z_vec)

SUITE["Basis Creation"] = @benchmarkable monomials(x, 3)
end

BenchMonomials.SUITE
