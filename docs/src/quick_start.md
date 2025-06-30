# [Quick Start](@id quick-start)

The non-commutative Broyden banded function is a generalization of the classical Broyden banded function to non-commuting variables. It is often used in optimization and numerical analysis to test the performance of algorithms. We will use it as an example.

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

The variables $x_i$ are non-commuting. You may think of them as matrices or operators that is assigned with a representation. It's possible puts constraints on the variables in the form of polynomials of equalities and inequalities. For example, we may require

```math
1 - x_i^2 \geq 0 \quad \text{and} \quad x_i - \frac{1}{3} \geq 0 \quad \forall i \in [1,n]
```

Firstly, use [`NCTSSoS.PolyOpt`](@ref) object to represent this problem. Since the constraints are inequalities, we pass it to `ineq_constraints` argument to constructor [`NCTSSoS.polyopt`](@ref).

```@repl quick-start
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

```@repl quick-start

solver_config = SolverConfig(optimizer=Mosek.Optimizer; mom_order=3, cs_algo=MF(), ts_algo=AsIsElimination())

result = cs_nctssos(pop, solver_config)
@assert result.objective ≈ 0.0 atol=1e-5
```

This function is a good choice for benchmarking and applying term and correlative sparsity for several reasons:

1.  **Structured Sparsity**: The function exhibits a clear "banded" structure. Each component $f_i$ only depends on $x_{i-1}, x_i,$ and $x_{i+1}$. This inherent sparsity in variable dependencies makes it an ideal candidate for correlative sparsity, where we exploit the fact that many pairs of variables do not appear together in terms.
2.  **Scalability**: The problem size (number of variables $n$) can be easily scaled, allowing for testing the performance of sparsity techniques on problems of varying dimensions.
3.  **Controlled Non-commutativity**: While the variables are non-commuting, the interactions are localized. This allows for a focused study of how non-commutativity interacts with sparsity.
4.  **Polynomial Structure**: The function is a sum of squares of polynomials. When we expand $f$, the resulting polynomial will have a certain number of terms. Term sparsity techniques can be applied by identifying and utilizing the fact that not all possible monomials (up to a certain degree) will be present in this polynomial representation. For example, a term like $x_1 x_5$ would not appear if $n$ is large enough and the band is narrow, which can be exploited by term sparsity.

By applying correlative sparsity, we can decompose the problem based on the limited interaction between variables (e.g., $x_i$ only interacts with its immediate neighbors). By applying term sparsity, we can reduce the size of the optimization problem by only considering the monomials that actually appear in the polynomial $f$.


As expected, term sparsity and correlative sparsity techniques significantly reduced the time taken to optimize the non-commutative Broyden banded function while maintaining the accuracy of the solution. The results matched with  [wangExploitingTermSparsity2021](@cite).

To summarize, the workflow for solving polynomial optimization  