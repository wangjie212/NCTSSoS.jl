module BenchPolyOpt
using BenchmarkTools
using NCTSSoS, MosekTools

const SUITE = BenchmarkGroup()

function make_poly(n)
    f = 0.0
    for i = 1:n
        jset = max(1, i - 5):min(n, i + 1)
        jset = setdiff(jset, i)
        f += (2x[i] + 5 * x[i]^3 + 1)^2
        f -= sum([
            4x[i] * x[j] +
            10x[i]^3 * x[j] +
            2x[j] +
            4x[i] * x[j]^2 +
            10x[i]^3 * x[j]^2 +
            2x[j]^2 for j in jset
        ])
        f += sum([
            x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset
        ])
    end
	return f
end

n = 20 
@ncpolyvar x[1:n]
f = make_poly(n)
order = 3
cons = vcat([(1 - x[i]^2) for i = 1:n], [(x[i] - 1 / 3) for i = 1:n])

SUITE["Example1"] = @benchmarkable result = cs_nctssos(pop, solver_config) setup = (
    begin
        pop = polyopt(f; ineq_constraints=cons)
        solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=order,
            cs_algo=MF(), ts_algo=MMD())
    end
)
end

BenchPolyOpt.SUITE
