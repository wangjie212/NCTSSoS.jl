# `FastPolynomials.jl` Benchmark

Temporarily, performance benchmarks for `FastPolynomials.jl` vs `DynamicPolynomials.jl` are placed here.


## Variables related

### Variable Creation 
```julia
using BenchmarkTools
# using DynamicPolynomials
# using NCTSSoS.FastPolynomials

@btime @ncpolyvar x[1:10000]

# FastPolynomials.jl
968.500 μs (100003 allocations: 2.90 MiB)

# DynamicPolynomials.jl
11.823 ms (120003 allocations: 3.44 MiB)
```

### Monomials of degree `d` creation

```julia
using BenchmarkTools
using NCTSSoS.FastPolynomials: monomials
using DynamicPolynomials

@ncpolyvar x[1:10]

@btime monomials($x,3)

# FastPolynomials.jl
109.292 μs (15993 allocations: 826.41 KiB)
# DynamicPolynomials.jl
73.917 μs (4023 allocations: 820.97 KiB)
``` 

## Monomials Related

```julia
using BenchmarkTools
using DynamicPolynomials
using NCTSSoS.FastPolynomials

@ncpolyvar x[1:10]
var_vec = [x[1],x[2],x[2],x[1],x[3]]
z_vec = [10,20,2,0,3]

@btime Monomial($var_vec, $z_vec)

# FastPolynomials.jl
73.484 ns (10 allocations: 592 bytes)
# DynamicPolynomials.jl
1.708 ns (0 allocations: 0 bytes)

# Pruning of zeros in exponents caused slow down
```

## Polynomials Related

```julia
using BenchmarkTools
using DynamicPolynomials
using NCTSSoS.FastPolynomials

n = 5
@ncpolyvar x[1:n]
function poly_create(x,n)
f = 0.0
for i in 1:n
	jset = max(1, i - 5):min(n, i + 1)
	jset = setdiff(jset, i)
    f += (2x[i] + 5 * x[i]^3 + 1)^2
    f -= sum([
		4x[i] * x[j] +
		10x[i]^3 * x[j] +
		2x[j] +
		4x[i] * x[j]^2 +
		10x[i]^3 * x[j]^2 +
		2x[j]^2 for j in jset])
		f += sum([
            x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset
        ])
end        
end

@btime poly_create($x, $n)
# FastPolynomials.jl
929.708 μs (29893 allocations: 1.81 MiB)

# DynamicPolynomials.jl
307.125 μs (18291 allocations: 939.89 KiB)
```

```julia
using BenchmarkTools
using DynamicPolynomials
using NCTSSoS.FastPolynomials

@ncpolyvar x[1:5]

p1 = x[1]^2 + 2x[2]*x[3] + 3x[3]*x[1] + 4x[4]^100*x[5]^2

@btime variables($p1)

# FastPolynomials.jl
524.869 ns (34 allocations: 3.28 KiB)
# DynamicPolynomials.jl
1.500 ns (0 allocations: 0 bytes)
```

## Comparison

### `in` operator

```julia
using BenchmarkTools
using NCTSSoS.FastPolynomials
using DynamicPolynomials

@ncpolyvar x[1:10000]

vars_vec = rand(x,5000)

@btime $x[1] in $vars_vec

# FastPolynomials.jl
26.541 μs (0 allocations: 0 bytes)

# DynamicPolynomials.jl
550.358 ns (0 allocations: 0 bytes)
```

## Real Problem

```julia 
using BenchmarkTools
using MosekTools, NCTSSoS
using Profile 

n = 5
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
    push!(cons, 1 - x[i]^2)
    push!(cons, x[i] - 1 / 3)
end
pop =  PolyOpt(f; constraints = cons);
solver_config = SolverConfig(optimizer=Mosek.Optimizer; mom_order=3, cs_algo=MF(), ts_algo=MMD())

@btime cs_nctssos($pop, $solver_config);

# NCTSSoS.jl
  8.963 s (357949474 allocations: 15.81 GiB)
# NCTSSOS.jl
  1.059 s (11859805 allocations: 478.06 MiB)

Profile.clear()
@profile result_cs_ts = cs_nctssos(pop, solver_config)
Profile.print(mincount=100, format=:flat)
```

## State Polynomial Problem
```julia
using BenchmarkTools
using MosekTools, NCTSSoS
using Profile 
using NCTSSoS.FastPolynomials: ς

@ncpolyvar x[1:3] y[1:3]
cov(a, b) = 1.0 * ς(x[a] * y[b]) - 1.0 * ς(x[a]) * ς(y[b])
sp = cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) + cov(3, 1) - cov(3, 2)

spop = StatePolyOpt(sp; is_unipotent=true, comm_gps=[x[1:3], y[1:3]])

solver_config = SolverConfig(; optimizer=Mosek.Optimizer, mom_order=2)

@btime cs_nctssos($spop, $solver_config)

Profile.clear()
@profile for _ in 1:20 result_cs_ts = cs_nctssos(spop, solver_config) end
Profile.print(mincount=100, format=:tree)
```


