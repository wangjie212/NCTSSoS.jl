using BenchmarkTools

SUITE = BenchmarkGroup()


for file in readdir(@__DIR__)
    # if startswith(file, "bench_") && endswith(file, ".jl")
    if startswith(file, "bench_") && endswith(file, "pop.jl")
        SUITE[file[(length("bench_")+1):(end-length(".jl"))]] = include(file)
    end
end
