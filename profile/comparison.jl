using BenchmarkTools

using DynamicPolynomials
@ncpolyvar x y z

mono1 = x^3*y
mono2 = x^2*y^2
mono3 = x^2*z*y


#   27.886 ns (0 allocations: 0 bytes)
#   50.109 ns (0 allocations: 0 bytes)
@btime cmp(mono1, mono2)

@which cmp(mono1, mono2)

#   377.702 ns (24 allocations: 1.12 KiB)
#   423.578 ns (24 allocations: 1.12 KiB)
@btime cmp(mono1, mono3)

using BenchmarkTools
using NCTSSoS.FastPolynomials

@ncpolyvar x y z

mono1 = Monomial([x, y], [3, 1])  
mono2 = Monomial([x, y], [2, 2])
mono3 = Monomial([x, z, y], [2, 1, 1])

#   51.244 ns (2 allocations: 64 bytes) mac mini before opt
#   48.878 ns (0 allocations: 0 bytes) mac pro
@btime cmp(mono1, mono2)

#   104.322 ns (2 allocations: 64 bytes)
@btime cmp(mono1, mono3)

@ncpolyvar x[1:100000]

var_coll = sort([x[i] for i in 1:99999])

# BenchmarkTools.Trial: 10000 samples with 960 evaluations per sample.
#  Range (min … max):  86.849 ns …  1.148 μs  ┊ GC (min … max): 0.00% … 89.56%
#  Time  (median):     95.095 ns              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   96.888 ns ± 26.765 ns  ┊ GC (mean ± σ):  0.77% ±  2.62%

#           ▃  ▃█▅                                               
#   ▃▇▆▃▂▃▆▇██████▇▄▄▄▃▃▂▂▂▂▂▂▂▂▂▂▁▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
#   86.8 ns         Histogram: frequency by time         127 ns <

#  Memory estimate: 32 bytes, allocs estimate: 1.

@benchmark searchsortedfirst(var_coll, x[100000])	