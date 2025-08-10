using NCTSSoS, MosekTools, Test
# TODO: Export FastPolynomials.StatePolynomial
using NCTSSoS.FastPolynomials: StatePolynomial

n = 5

function pauli_operators(n::Int)
    @ncpolyvar x[1:n]  # Pauli x operators
    @ncpolyvar y[1:n]  # Pauli y operators
    @ncpolyvar z[1:n]  # Pauli z operators
    equalities = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:n])
    comm_gps = [[x[i], y[i], z[i]] for i in 1:n]
    return x, y, z, comm_gps, equalities
end

# extra constraints
cov(a, b) = 1.0 * ς(a * b) * one(Monomial)
function rydberg_blockade(z)
    n = length(z)
    return [cov(z[i], z[i+1]) for i in 1:n-1]
end

x, y, z, comm_gps, equalities = pauli_operators(n)
f = sum(z[i] * z[i+1] for i in 1:n-1)  # objective function

# TODO: cast to state poly automatically!!!!
f = NCStatePolynomial(f)
equalities = NCStatePolynomial.(equalities)
blockade = rydberg_blockade(z)

pop = polyopt(             # the optimization problem
        f;
        comm_gps,   # items in the different groups commute
        eq_constraints = [equalities; blockade],
        is_unipotent=true  # unipotent: X^2 = 1
    )

solver_config = SolverConfig(;
    optimizer=Mosek.Optimizer,  # the solver backend
    order=1                    # the order of the moment matrix
)
result = cs_nctssos(pop, solver_config)
result.objective  # the upper bound of the CHSH inequality