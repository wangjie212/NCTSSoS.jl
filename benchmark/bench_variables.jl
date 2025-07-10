module BenchVariables
using BenchmarkTools
using NCTSSoS.FastPolynomials

const SUITE = BenchmarkGroup()

N = 10000

SUITE["Variables Creation"] = @benchmarkable @ncpolyvar x[1:N]

@ncpolyvar x[1:N]

SUITE["Variable Test `in`"] = @benchmarkable Base.in(x[1], vars_vec) setup = vars_vec = sort(rand(x, 5000))
end

BenchVariables.SUITE
