# [Quick Start](@id quick-start)

The non-commutative Broyden banded function is a generalization of the classical
Broyden banded function to non-commuting variables. It is often used in
optimization and numerical analysis to test the performance of algorithms. We
will use it as an example.

The function is defined as:
```math
f(x_1, \dots, x_n) = \sum_{i=1}^n f_i(x_1, \dots, x_n)^2
```
where
```math
f_i(x_1, \dots, x_n) = 2x_i + 5x_i^3 + 1 - \sum_{j \in J_i} (x_j +x_j)^2
```
with
```math
J_i = \{j | j \neq i, \max(1, i-5) \leq j \leq \min(n, i+1)\}
```

The variables $x_i$ are non-commuting. You may think of them as matrices or
operators that is assigned with a representation. It's possible puts constraints
on the variables in the form of polynomials of equalities and inequalities. For
example, we may require

```math
1 - x_i^2 \geq 0 \quad \text{and} \quad x_i - \frac{1}{3} \geq 0 \quad \forall i \in [1,n]
```

Firstly, use [`NCTSSoS.PolyOpt`](@ref) object to represent this problem. Since
the constraints are inequalities, we pass it to `ineq_constraints` argument to
constructor [`NCTSSoS.polyopt`](@ref).

```julia quick-start
using NCTSSoS, MosekTools

function broyden_banded(n::Int)
	@ncpolyvar x[1:n]
	f = 0.
	for i = 1:n
	    jset = max(1, i-5) : min(n, i+1)
	    jset = setdiff(jset, i)
	    g = sum(x[j] + x[j]^2 for j in jset)
	    f += (2*x[i] + 5*x[i]^3 + 1 - g)^2
	end
	ineq_cons = [[1 - x[i]^2 for i in 1:n];[x[i] - 1/3 for i in 1:n]]

	return polyopt(f; ineq_constraints=ineq_cons)
end

pop = broyden_banded(6)
```
[Polynomial Optimization](@ref polynomial-optimization) is NP-hard, therefore it
is considered impossible to solve them efficiently in general. However, it is
possible to relax the problem into [Semidefinite Programming](@ref
semidefinite-programming). The solution of the Semidefinite program will be the
lower/upper bound of the original minimization/maximization polynomial
optimziation problem.

The relaxation can be done in two different forms, [moment relaxation](@ref
moment-problem) and [sum of hermitian square relaxation](@ref sohs-problem). The
relaxation is tight in the limit of a parameter, moment order, reaching
infinity. However, you may get lucky and be able to sovle the problem at finite
moment order.

Information about the relaxation is encoded in [`NCTSSoS.SolverConfig`](@ref).
Besides moment order, it also needs to be provided with a [SDP Solver](@ref
overview-of-optimizers). Different optimizers may have different performance
characteristics and may be more or less suitable for a given problem.

```julia quick-start
solver_config = SolverConfig(optimizer=Mosek.Optimizer;mom_order=3)
```

Finally, we are ready to solve this problem. This is accomplishded with [`NCTSSoS.cs_nctssos`](@ref).

```julia quick-start
@time result = cs_nctssos(pop, solver_config)

@assert result.objective ≈ 0.0 atol=1e-5 # find correct value
```

Although we have reached a tight bound, time to solve the problem can still be
significant. As a remedy, [Sparsities](@ref sparsities) can be utilized to
reduce the problem size. This is achieved by supplying an
[`EliminationAlgorithm`](https://algebraicjulia.github.io/CliqueTrees.jl/stable/api/#Elimination-Algorithms)
to [`NCTSSoS.SolverConfig`](@ref).

```julia quick-start
solver_config = SolverConfig(optimizer=Mosek.Optimizer; mom_order=3, cs_algo=MF(), ts_algo=AsIsElimination())

@time result = cs_nctssos(pop, solver_config)

@assert result.objective ≈ 0.0 atol=1e-5 # ideally, this has to fail
```

As expected, the time taken to solve the problem has been significantly reduced.
However, the result is no longer a tight lower bound. Sparsity is itself a kind
of relaxation. Luckily, we may tighten the relaxation in the [Term
Sparsity](@ref term-sparsity) sense, this is done with
[`NCTSSoS.cs_nctssos_higher`](@ref).

```julia quick-start

@time result_higher = cs_nctssos_higher(pop, solver_config)

@assert result_higher.objective ≈ 0.0 atol=1e-5 # ideally, this should approach the optimal value
```

## Workflow

To summarize, the workflow for solving polynomial optimization can be summarized as

![`Workflow for solving Polynomial Optimization problem`](assets/overall_workflow.typ.svg)

If you would like to understand more please refer to [examples section] in the document.

## Why is Broyden Banded Function Desirable

This function is a good choice for benchmarking and applying term and correlative sparsity for several reasons:

1.  **Structured Sparsity**: The function exhibits a clear "banded" structure. Each component $f_i$ only depends on $x_{i-1}, x_i,$ and $x_{i+1}$. This inherent sparsity in variable dependencies makes it an ideal candidate for correlative sparsity, where we exploit the fact that many pairs of variables do not appear together in terms.
2.  **Scalability**: The problem size (number of variables $n$) can be easily scaled, allowing for testing the performance of sparsity techniques on problems of varying dimensions.
3.  **Controlled Non-commutativity**: While the variables are non-commuting, the interactions are localized. This allows for a focused study of how non-commutativity interacts with sparsity.
4.  **Polynomial Structure**: The function is a sum of squares of polynomials. When we expand $f$, the resulting polynomial will have a certain number of terms. Term sparsity techniques can be applied by identifying and utilizing the fact that not all possible monomials (up to a certain degree) will be present in this polynomial representation. For example, a term like $x_1 x_5$ would not appear if $n$ is large enough and the band is narrow, which can be exploited by term sparsity.

By applying correlative sparsity, we can decompose the problem based on the limited interaction between variables (e.g., $x_i$ only interacts with its immediate neighbors). By applying term sparsity, we can reduce the size of the optimization problem by only considering the monomials that actually appear in the polynomial $f$.

As expected, term sparsity and correlative sparsity techniques significantly reduced the time taken to optimize the non-commutative Broyden banded function while maintaining the accuracy of the solution. The results matched with  [wangExploitingTermSparsity2021](@cite).
