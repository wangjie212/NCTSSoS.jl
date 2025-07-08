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

n = 12
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

# current record
# |                                      | main             | dirty             | main / dirty   |
# |:-------------------------------------|:----------------:|:-----------------:|:--------------:|
# | monomials/Basis Creation             | 0.105 ± 0.004 ms | 0.105 ± 0.0039 ms | 0.999 ± 0.053  |
# | monomials/Compare different degree   | 0.001 ± 41 ns    | 0.001 ± 41 ns     | 1 ± 5.8e+04    |
# | monomials/Compare same degree        | 0.001 ± 41 ns    | 0.001 ± 41 ns     | 1 ± 5.8e+04    |
# | monomials/Monomial Creation          | 0.084 ± 0.042 μs | 0.084 ± 0.042 μs  | 1 ± 0.71       |
# | polynomials/Polynomial Creation      | 0.806 ± 0.023 ms | 0.823 ± 0.016 ms  | 0.98 ± 0.034   |
# | polynomials/Polynomial get variables | 0.542 ± 0.043 μs | 0.542 ± 0.083 μs  | 1 ± 0.17       |
# | pop/Example1                         | 14.2 s           | 3.46 ± 0.078 s    | 4.1            |
# | variables/Variable Test `in`         | 23.8 ± 2 μs      | 0.166 ± 0.042 μs  | 143 ± 38       |
# | variables/Variables Creation         | 5.44 ± 0.78 ms   | 5.45 ± 0.89 ms    | 0.998 ± 0.22   |
# | time_to_load                         | 0.915 ± 0.0065 s | 1.11 ± 0.0019 s   | 0.827 ± 0.0061 |