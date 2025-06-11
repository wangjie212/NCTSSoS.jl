# Correlative and Term Sparsity

The non-commutative Broyden banded function is a generalization of the classical Broyden banded function to non-commuting variables. It is often used in optimization and numerical analysis to test the performance of algorithms.

The function is defined as:
$$f(x_1, \dots, x_n) = \sum_{i=1}^n f_i(x_1, \dots, x_n)^2$$
where
$$f_i(x_1, \dots, x_n) = 2x_i + 5x_i^3 + 1 - \sum_{j \in J_i} (X_j +X_j)^2$$
with $J_i = \{j | j \neq i, max(1, i-5) \leq j \leq min(n, i+1)\}$. The variables $x_i$ are non-commuting.

This function is a good choice for benchmarking and applying term and correlative sparsity for several reasons:

1.  **Structured Sparsity**: The function exhibits a clear "banded" structure. Each component $f_i$ only depends on $x_{i-1}, x_i,$ and $x_{i+1}$. This inherent sparsity in variable dependencies makes it an ideal candidate for correlative sparsity, where we exploit the fact that many pairs of variables do not appear together in terms.
2.  **Scalability**: The problem size (number of variables $n$) can be easily scaled, allowing for testing the performance of sparsity techniques on problems of varying dimensions.
3.  **Controlled Non-commutativity**: While the variables are non-commuting, the interactions are localized. This allows for a focused study of how non-commutativity interacts with sparsity.
4.  **Polynomial Structure**: The function is a sum of squares of polynomials. When we expand $f$, the resulting polynomial will have a certain number of terms. Term sparsity techniques can be applied by identifying and utilizing the fact that not all possible monomials (up to a certain degree) will be present in this polynomial representation. For example, a term like $x_1 x_5$ would not appear if $n$ is large enough and the band is narrow, which can be exploited by term sparsity.

By applying correlative sparsity, we can decompose the problem based on the limited interaction between variables (e.g., $x_i$ only interacts with its immediate neighbors). By applying term sparsity, we can reduce the size of the optimization problem by only considering the monomials that actually appear in the polynomial $f$.

## Non-commutative Broyden Banded Function

```julia 
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

	pop = PolyOpt(f)

	solver_config = SolverConfig(optimizer=Mosek.Optimizer; mom_order=3, cs_algo=MF(), ts_algo=AsIsElimination())

	result = cs_nctssos(pop, solver_config)
	@assert result.objective ≈ 0.0 atol=1e-5
end

elapsed_time = [@elapsed broyden_banded(n) for n in[20,40,60,80,100,200,300,400,500,600,700,800,900,1000]]
@show elapsed_time
```

```julia

```

As expected, term sparsity and correlative sparsity techniques significantly reduced the time taken to optimize the non-commutative Broyden banded function while maintaining the accuracy of the solution. The results matched with [^Wang].

[^Wang]: Wang, J. and Magron, V., 2021. Exploiting term sparsity in noncommutative polynomial optimization. Computational Optimization and Applications, 80(2), pp.483-521.