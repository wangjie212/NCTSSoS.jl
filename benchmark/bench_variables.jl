module BenchVariables
using BenchmarkTools
using NCTSSoS.FastPolynomials

const SUITE = BenchmarkGroup()


N = 10000
@ncpolyvar x[1:N]

SUITE["Variables Creation"] = @benchmarkable @ncpolyvar x[1:N]

vars_vec = rand(x, 5000)

SUITE["Variable Test `in`"] = x[1] in vars_vec

end

BenchVariables.SUITE
