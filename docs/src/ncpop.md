# Noncommutative polynomial optimization

Noncommutative polynomial optimization concerns minimizing the smallest eigenvalue or the trace of a noncommutative polynomial subject to a tuple of noncommutative polynomial inequality constraints and equality constraints, which in general takes the form:

$$\mathrm{inf}_{\mathbf{x}\in\mathcal{B}(\mathcal{H})^n}\ \lambda_{\min}(f(\mathbf{x}))\ (\text{or } \mathrm{tr}(f(\mathbf{x}))) \ \text{ s.t. }\ g_1(\mathbf{x})\ge0,\ldots,g_m(\mathbf{x})\ge0,h_1(\mathbf{x})=0,\ldots,h_{\ell}(\mathbf{x})=0,$$

where $f,g_1,\ldots,g_m,h_1,\ldots,h_{\ell}\in\mathbb{R}\langle\mathbf{x}\rangle$ are noncommutative polynomials in noncommuting variables $\mathbf{x}$, and $\mathcal{H}$ is an (infinite dimensional) seperable Hilbert space.

To illustrate how to solve a noncommutative polynomial optimization problem with NCTSSOS, let us consider the following simple example.

```Julia
using NCTSSoS
using MosekTools
@ncpolyvar x[1:2]
f = 2.0 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
ineq = [4.0 - x[1]^2 - x[2]^2]
eq = [x[1]*x[2] + x[2]*x[1] - 2.0]
pop = PolyOpt(f; constraints=[eq;ineq], is_equality=[true,false])
d = 2 # set the relaxation order
result = cs_nctssos(pop, SolverConfig(optimizer=Mosek.Optimizer;mom_order=2)) # compute the first TS step of the NCTSSOS hierarchy
result_higher = cs_nctssos_higher(pop,result,SolverConfig(optimizer=Mosek.Optimizer; mom_order=2)) # compute higher TS steps of the NCTSSOS hierarchy
```

## Correlative sparsity
The following is an example where one exploits correlative sparsity and term sparsity simultaneously.

```Julia
using NCTSSoS
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

cons = vcat([(1 - x[i]^2) for i in 1:n], [(x[i] - 1 / 3) for i in 1:n])

pop = PolyOpt(f; constraints=cons)
d = 3 # set the relaxation order
result = cs_nctssos(pop, SolverConfig(optimizer=Mosek.Optimizer; mom_order=d, cs_algo=MF(), ts_algo=MMD())) # compute the first TS step of the CS-NCTSSOS hierarchy
result_higher = cs_nctssos_higher(pop,result,SolverConfig(optimizer=Mosek.Optimizer; mom_order=d, cs_algo=MF(), ts_algo=MMD())) # compute higher TS steps of the CS-NCTSSOS hierarchy
```

## Methods
```@docs
cs_nctssos
cs_nctssos_higher
```

### References

1. [Exploiting Term Sparsity in Noncommutative Polynomial Optimization](https://arxiv.org/abs/2010.06956), 2021.    
2. [Sparse polynomial optimization: theory and practice](https://arxiv.org/abs/2208.11158), 2023.   
3. [Optimization of polynomials in non-commuting variables](https://link.springer.com/content/pdf/10.1007/978-3-319-33338-0.pdf), 2016. 