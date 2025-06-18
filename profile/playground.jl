using DynamicPolynomials
using DynamicPolynomials: Monomial
using BenchmarkTools, Profile

@ncpolyvar x[1:10]

exponents = [rand(0:4, 10) for _ in 1:1000]

#   184.542 μs (6494 allocations: 539.38 KiB)
@btime monos = [Monomial(x, deepcopy(exponents[j])) for j in 1:1000]

#   142.875 μs (5998 allocations: 345.89 KiB)
@btime mono_vec = MonomialVector(x, deepcopy(exponents))