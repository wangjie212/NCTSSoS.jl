module BenchVariables
	using BenchmarkTools
	using NCTSSoS.FastPolynomials

	const SUITE = BenchmarkGroup()

	SUITE["Variables Creation"] = @benchmarkable @ncpolyvar x[1:10000]
end

BenchVariables.SUITE
