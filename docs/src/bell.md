# Bell inequalities

## Linear bell inequalities

The CHSH inequality:

```Julia
using NCTSSoS
@ncpolyvar x[1:2]
@ncpolyvar y[1:2]
f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
pop = PolyOpt(f; comm_gp= Set(x), is_unipotent=true)

solver_config = SolverConfig(optimizer=Clarabel.Optimizer; mom_order=1)

result = cs_nctssos(pop, solver_config)
```

The I_3322 inequality:

```Julia
using NCTSSoS
@ncpolyvar x[1:3]
@ncpolyvar y[1:3]
f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) + x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]
pop = PolyOpt(-f; comm_gp= Set(x), is_projective=true)

solver_config = SolverConfig(optimizer=Clarabel.Optimizer; mom_order=3)

result = cs_nctssos(pop, solver_config)
```