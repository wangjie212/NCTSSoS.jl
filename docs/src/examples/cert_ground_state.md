# Certifying Ground State

Finding the ground state of a quantum system is a fundamental problem in quantum
mechanics [wang2024Certifying](@cite). Variational methods are commonly used to
approximate the ground state. Due to the variational nature of these methods,
only an upper bound can be obtained [kull2024Lower](@cite). Polynomial
optimization techniques provides a way to find the lower bound of the ground
state energy. With such lower bound, estimates of properties of the ground
state, correlations functions, structure factors and order parameters, can also
be obtained. We provide examples on 1D and 2D Heisenberg models.

In general, we consider the following Hamiltonian:
```math
H = \frac{1}{4} \sum_{i < j} J_{ij} \sum_{ a \in \{x,y,z\}} \sigma_i^a \sigma_j^a
```

## 1D Heisenberg Model with Nearest Neighbor Interaction

Firstly, let's consider the simplest case of 1D Heisenberg chain with nearest
neighbor interaction and periodic boundary condition.

```julia 1D-Heisenberg
using NCTSSoS, MosekTools
N = 6
@ncpolyvar x[1:N] y[1:N] z[1:N]

ham = sum(ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [x, y, z] for i in 1:N)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

pop = polyopt(
            ham;                                        # the Hamiltonian
            eq_constraints=eq_cons,                     # anti-commutation relation between Pauli Operators
            comm_gps=[[x[i], y[i], z[i]] for i in 1:N], # commutation relation between Pauli Operators
            is_unipotent=true                           # Pauli operators square to identity
            )

solver_config = SolverConfig(
                    optimizer=Mosek.Optimizer,          # the solver backend
                    order=2,                        # moment matrix order
                    ts_algo = MMD(),                    # term sparsity algorithm
                    )

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(
            pop,                                        # Polynomial Optimization Problem
            res,                                        # Solution of First Order Term Sparsity Iteration
            solver_config                               # Solver Configuration
        )
res.objective / N
```

```julia
-0.46712927163730084
```

The returned result matches the actual ground state energy $-0.467129$ to $6$
digits. [wang2024Certifying](@cite)


## 1D Heisenberg Model with next nearest neighbor interaction

Polynomial Optimization framework is quite general. Almost no modification is
required to handle more complex Hamiltonian. 1D Heisenberg Model with geometric
frustration induced by next nearest neighbor interaction can be solved as:

```julia geom-frustration
using NCTSSoS, MosekTools
N = 6
J1 = 1.0                            # Nearest Neighbor Interaction
J2 = 0.2                            # Next Nearest Neighbor Interaction
@ncpolyvar x[1:N] y[1:N] z[1:N]

ham = sum(ComplexF64(J1 / 4) * op[i] * op[mod1(i + 1, N)] + ComplexF64(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

pop = polyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=2, ts_algo = MMD())

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(pop,res,solver_config)
res.objective/N
```
```julia
-0.4270083225088391
```

We are able to obtain the ground state energy of $-0.4270083225302217$, accurate
to $6$ digits!

## 2D Square Lattice

Extending Heisenberg model to $2$-D case is also straightforward.

```julia
using NCTSSoS, MosekTools
Nx = 3
Ny = 3
N = Nx * Ny
J1 = 1.0
J2 = 0.0
@ncpolyvar x[1:N] y[1:N] z[1:N]

LI = LinearIndices((1:Nx, 1:Ny))

ham = sum(ComplexF64(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(i, mod1(j + 1, Ny))]] + ComplexF64(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), j)]] + ComplexF64(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j + 1, Ny))]] + ComplexF64(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j - 1, Ny))]] for op in [x, y, z] for i in 1:Nx for j in 1:Ny)

eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

pop = polyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=3, cs_algo=MF(), ts_algo=MMD())

res = cs_nctssos(pop, solver_config)

res = cs_nctssos_higher(pop,res,solver_config)
res = cs_nctssos_higher(pop,res,solver_config)

res.objective / N
```

## 1D Heisenberg Model with Next Nearest Neighbor Interaction

Now that we can find the lower bound of energy, let's try to bound the operator
value with respect to the ground state.

To do so, we add two inequality constraints and change the objective value to
the target operator `x[1]x[2]`

```julia
using NCTSSoS,MosekTools

SOLVER = Mosek.Optimizer

energy_upper_bounds = [-0.44667162694019086,
    -0.42700832253022214,
    -0.40833333333333327,
    -0.390893734117895,
    -0.375,
    -0.4000000000000003,
    -0.42500000000000004,
    -0.45,
    -0.4750000000000003,
    -0.5000000000000001,
    -0.5249999999999998,
    -0.55,
    -0.5749999999999996,
    -0.6000000000000002,
    -0.6249999999999999,
    -0.6500000000000002,
    -0.6750000000000003,
    -0.7000000000000003,
    -0.7249999999999998,
    -0.7500000000000003]


s0s1_vals = [-0.15558685225793648,
    -0.1551397955033638,
    -0.1542145593869776,
    -0.1525967357708984,
    -0.24998818885551125,
    -0.08333333333445211,
    -0.08333333333333304,
    -0.08333333333333384,
    -0.08333333333333358,
    -0.0833333333333331,
    -0.08333333333333363,
    -0.08333333333333394,
    -0.08333333333333363,
    -0.08333333333333329,
    -0.08333333333333329,
    -0.08333333333333329,
    -0.08333333333333227,
    -0.08333333333333347,
    -0.08333333333333348,
    -0.0833333333333339]

function main(N::Int, J2s, energy_upper_bounds, s0s1_vals)
    J1 = 1.0

    op_upper_bounds = Float64[]
    op_lower_bounds = Float64[]
    energy_lower_bounds = Float64[]

    for (idx, J2) in enumerate(J2s)
        @ncpolyvar x[1:N] y[1:N] z[1:N]

        ham = sum(ComplexF64(J1 / 4) * op[i] * op[mod1(i + 1, N)] + ComplexF64(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

        eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

        pop = PolyOpt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

        solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())

        res = cs_nctssos(pop, solver_config)

        res = cs_nctssos_higher(pop, res, solver_config)

        energy_lower_bound = res.objective

        push!(energy_lower_bounds, energy_lower_bound)

        @ncpolyvar x[1:N] y[1:N] z[1:N]

        obj = one(ComplexF64) * x[1] * x[2]

        ham = sum(ComplexF64(J1 / 4) * op[i] * op[mod1(i + 1, N)] + ComplexF64(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

        eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

        ineq_cons = -1 .* [energy_upper_bounds[idx]*N - ham, ham - energy_lower_bound]

        pop = PolyOpt(obj; eq_constraints=eq_cons, ineq_constraints=ineq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

        solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())

        res = cs_nctssos(pop, solver_config)

        res = cs_nctssos_higher(pop, res, solver_config)

        operator_lower_bound = res.objective

        push!(op_lower_bounds, operator_lower_bound)

        pop = PolyOpt(-obj; eq_constraints=eq_cons, ineq_constraints=ineq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

        solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())

        res = cs_nctssos(pop, solver_config)

        res = cs_nctssos_higher(pop, res, solver_config)

        operator_upper_bound = -res.objective
        push!(op_upper_bounds, operator_upper_bound)
    end

    open("N$(N)J1J2_pop.txt", "w") do file
        writedlm(file, hcat(op_upper_bounds, op_lower_bounds, energy_upper_bounds, energy_lower_bounds, s0s1_vals), ' ')
    end
end

main(6, 0.1:0.1:2.0, energy_upper_bounds, s0s1_vals)
```

We produce the following result, which does not align with [wang2024Certifying](@cite).
