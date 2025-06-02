using BenchmarkTools

using DynamicPolynomials
@ncpolyvar x y z

mono1 = x^3*y
mono2 = x^2*y^2
mono3 = x^2*z*y


#   27.886 ns (0 allocations: 0 bytes)
@btime cmp(mono1, mono2)
#   377.702 ns (24 allocations: 1.12 KiB)
@btime cmp(mono1, mono3)

using BenchmarkTools
using NCTSSoS.FastPolynomials

@ncpolyvar x y z

mono1 = Monomial([x, y], [3, 1])  
mono2 = Monomial([x, y], [2, 2])
mono3 = Monomial([x, z, y], [2, 1, 1])

#   51.244 ns (2 allocations: 64 bytes)
@btime cmp(mono1, mono2)
@btime cmp(mono1, mono3)