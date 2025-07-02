# Tracial Polynomial Optimization


```julia
using NCTSSoS, MosekTools
using NCTSSoS.FastPolynomials:tr, Monomial
@ncpolyvar x[1:3]

p = (tr(x[1] * x[2] * x[3]) + tr(x[1] * x[2]) * tr(x[3])) * one(Monomial)

spop = polyopt(p; is_projective=true, comm_gps=[x])

solver_config = SolverConfig(; optimizer=SOLVER, mom_order=2)

result = cs_nctssos(spop, solver_config)

@test result.objective ≈ -0.046717378455438933 atol = 1e-6

solver_config = SolverConfig(; optimizer=Mosek.Optimizer, mom_order=3)

result = cs_nctssos(spop, solver_config)

@test result.objective ≈ -0.03124998978001017 atol = 1e-6
```

```julia
using NCTSSoS, MosekTools
using NCTSSoS.FastPolynomials:tr, Monomial

	@ncpolyvar x[1:2] y[1:2]

    p = -1.0 * tr(x[1] * y[1]) - 1.0 * tr(x[1] * y[2]) - 1.0 * tr(x[2] * y[1]) + 1.0 * tr(x[2] * y[2])

	tpop = polyopt(p * one(Monomial); is_unipotent=true)

	solver_config = SolverConfig(; optimizer=Mosek.Optimizer, mom_order=1, ts_algo=MaximalElimination())

	result = cs_nctssos(tpop, solver_config)

    @test result.objective ≈ -2.8284271157283083 atol = 1e-5
````
