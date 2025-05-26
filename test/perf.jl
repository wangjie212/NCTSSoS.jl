using BenchmarkTools
using MosekTools, NCTSSoS

@testset "Bryoden Banded Polynomial" begin
    n = 10
    @ncpolyvar x[1:n]
    f = 0.0
    for i = 1:n
        jset = max(1, i - 5):min(n, i + 1)
        jset = setdiff(jset, i)
        f += (2x[i] + 5 * x[i]^3 + 1)^2
        f -= sum([4x[i] * x[j] + 10x[i]^3 * x[j] + 2x[j] + 4x[i] * x[j]^2 + 10x[i]^3 * x[j]^2 + 2x[j]^2 for j in jset])
        f += sum([x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset])
    end

	cons = typeof(f)[]
    for i = 1:n
        push!(pop, 1 - x[i]^2)
        push!(pop, x[i] - 1 / 3)
    end

    pop =  PolyOpt(f);

    result_cs_ts = cs_nctssos(pop, SolverConfig(optimizer=Mosek.Optimizer; cs_algo=MF(), ts_algo=MMD()))

    result_cs_ts_higher = cs_nctssos_higher(pop, result_cs_ts, SolverConfig(optimizer=Clarabel.Optimizer; cs_algo=MF(), ts_algo=MMD()))
end
