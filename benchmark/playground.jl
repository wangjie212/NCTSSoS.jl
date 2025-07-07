using BenchmarkTools
using NCTSSoS.FastPolynomials 
using NCTSSoS.FastPolynomials: monomial, monomials

@ncpolyvar x[1:10]

@benchmark degree(a) setup = (a = monomial([x[1], x[2]], [2, 1]))

@benchmark cmp(a, b) setup = (a = monomial([x[1], x[2]], [2, 1]); b = monomial([x[1], x[2]], [1, 1]))

@benchmark cmp(a, b) setup = (a = monomial([x[1], x[2], x[1]], [2, 1, 1]); b = monomial([x[1], x[2]], [3, 1]))