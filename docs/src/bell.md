# Bell inequalities

Bell inequalities are mathematical expressions that test whether the predictions of quantum mechanics can be explained by local hidden variable theories. They were first introduced by John Stewart Bell in 1964 and have since become fundamental tools in quantum information theory and quantum foundations.
A Bell inequality is typically expressed as a linear combination of expectation values of observables, with bounds that differ between classical and quantum theories. In the classical case, these inequalities must be satisfied if the system can be described by local hidden variables. However, quantum mechanics can violate these inequalities, demonstrating the non-local nature of quantum correlations.
The general form of a Bell inequality can be written as:

$$\sum_{i,j} c_{ij} \langle A_i B_j \rangle \leq C$$

where $A_i$ and $B_j$ are observables measured by two parties (traditionally called Alice and Bob), $c_{ij}$ are real coefficients, and $C$ is the classical bound. Quantum mechanics can violate this inequality, with the maximum violation known as the quantum bound.

## Linear Bell inequalities

### CHSH inequality
The most famous Bell inequality is the CHSH (Clauser-Horne-Shimony-Holt) inequality, which involves two parties, each measuring two observables. For binary observables $A_1, A_2$ measured by Alice and $B_1, B_2$ measured by Bob. We define the objective function as:

$$f(A_1, A_2, B_1, B_2) = \langle A_1B_1 \rangle + \langle A_1B_2 \rangle + \langle A_2B_1 \rangle - \langle A_2B_2 \rangle$$

The CHSH inequality is then given by $$f(A_1, A_2, B_1, B_2) \leq 2$$, which must be satisfied by any local hidden variable theory. However, quantum mechanics can violate this inequality up to the value $2\sqrt{2}$, known as the Tsirelson bound. This violation demonstrates that quantum mechanics cannot be described by any local hidden variable theory.
The CHSH inequality is particularly important because it is the simplest non-trivial Bell inequality and has been experimentally verified numerous times, providing strong evidence for the non-local nature of quantum mechanics.

The upper bound of the CHSH inequality can be computed using the following code:

```@example chsh
using NCTSSoS, Clarabel

@ncpolyvar x[1:2]  # x[1] = A_1, x[2] = A_2
@ncpolyvar y[1:2]  # y[1] = B_1, y[2] = B_2
f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]  # objective function

pop = PolyOpt(             # the optimization problem
        f; 
        comm_gp= Set(x),   # commutative group
        is_unipotent=true  # the variables are unipotent
    )

solver_config = SolverConfig(;
    optimizer=Clarabel.Optimizer,  # the solver backend
    mom_order=1                    # the order of the moment matrix
)
result = cs_nctssos(pop, solver_config)
result.objective  # the upper bound of the CHSH inequality
```

Here, we first declare some operators as non-commutative variables, and then construct the optimization problem. In `PolyOpt` constructor,
- `comm_gp` argument specifies the commutative group of the variables, which means variables in different commutative groups commute with each other.
- `is_unipotent` argument specifies that the variables are unipotent, which means they are positive semidefinite.

Here, since the variables on different qubits commute with each other, we can group them into different commutative groups.

In the solver configuration, we use [Clarabel](https://github.com/oxfordcontrol/Clarabel.jl) as the semidefinite programming solver backend. It is an open-source solver for conic programs with quadratic objectives, and it uses the interior point method to solve the problem[^Goulart2024].

The resulting upper bound is very close to the theoretical exact value $2\sqrt{2} \approx 2.8284271247461903$ (accurate up to 7 decimals!!).

### $I_{3322}$ inequality

The $I_{3322}$ inequality is a more complex inequality that involves three parties, each measuring three observables. Let $A_1, A_2, A_3$ be the observables measured by Alice and $B_1, B_2, B_3$ be the observables measured by Bob. We define the objective function as[^Pal2010]:

```math
f(A_1, A_2, A_3, B_1, B_2, B_3) = \langle A_1(B_1+B_2+B_3) \rangle + \langle A_2(B_1+B_2-B_3) \rangle\\
+ \langle A_3(B_1-B_2) \rangle 
- \langle A_1 \rangle - 2\langle B_1 \rangle - \langle B_2 \rangle
```

In classical mechanics, the inequality $f(A_1, A_2, A_3, B_1, B_2, B_3) \leq 0$ must be satisfied. However, quantum mechanics can violate this inequality up to the value $0.25$, known as the Tsirelson bound. This violation demonstrates that quantum mechanics cannot be described by any local hidden variable theory.

The upper bound of the I3322 inequality can be computed using the following code:

```@example i3322
using NCTSSoS, Clarabel

@ncpolyvar x[1:3]
@ncpolyvar y[1:3]
f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) + x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]

pop = PolyOpt(-f; comm_gp= Set(x), is_projective=true)

solver_config = SolverConfig(optimizer=Clarabel.Optimizer; mom_order=2)

result = cs_nctssos(pop, solver_config)
result.objective
```

The resulting upper bound is close to the theoretical exact value $0.25$. By increasing the order of the moment matrix, this upper bound can be improved.

[^Goulart2024]: Goulart, P.J., Chen, Y., 2024. Clarabel: An interior-point solver for conic programs with quadratic objectives. https://doi.org/10.48550/arXiv.2405.12762
[^Pal2010]: Pál, K.F., Vértesi, T., 2010. Maximal violation of the I3322 inequality using infinite dimensional quantum systems. Phys. Rev. A 82, 022116. https://doi.org/10.1103/PhysRevA.82.022116
