using BenchmarkTools
using MosekTools, NCTSSoS
using Profile 

order = 3
n = 15 
@ncpolyvar x[1:n]
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

cons = vcat([(1 - x[i]^2) for i = 1:n], [(x[i] - 1 / 3) for i = 1:n])

pop = polyopt(f; ineq_constraints=cons)
solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=order,
    cs_algo=MF(), ts_algo=MMD())


@benchmark result = cs_nctssos($pop, $solver_config) 
# used to be 11 second
# BenchmarkTools.Trial: 2 samples with 1 evaluation per sample.
#  Range (min … max):  3.764 s …   3.801 s  ┊ GC (min … max): 11.23% … 11.84%
#  Time  (median):     3.783 s              ┊ GC (median):    11.54%
#  Time  (mean ± σ):   3.783 s ± 26.299 ms  ┊ GC (mean ± σ):  11.54% ±  0.43%

Profile.clear()
@profile result = cs_nctssos(pop, solver_config)
Profile.print(mincount=300)