abstract type OptimizationProblem end

@enum Objective EIGEN TRACE
# T: type of the coefficients
# OBJ: type of the objective function, one of `Objective`
struct PolyOpt{T,OBJ} <: OptimizationProblem
    objective::Polynomial{T}
    constraints::Vector{Polynomial{T}} # NOTE: assuming constraints are all simplified using comm_gp, is_unipotent, and is_projective
    is_equality::Vector{Bool} # which constraints are equality constraints
    variables::Vector{Variable}
    comm_gps::Vector{Set{Variable}} # Vectors of Set of variables that commutes with variables not in the same set
    is_unipotent::Bool # square to 1. Examples: Pauli Operators, SWAP Operators
    is_projective::Bool # X^2 = X. Is projective.
end

function PolyOpt(objective::Polynomial{T}; constraints=Any[], is_equality=fill(false, length(constraints)), comm_gps=Set{Variable}[], is_unipotent::Bool=false, is_projective::Bool=false, obj_type::Objective=EIGEN) where {T}
    @assert !(T <: Integer) "The polynomial coefficients can not be integers (not supported by JuMP solvers)."
    cons = collect(Polynomial{T}, constraints)
    is_eq = collect(Bool, is_equality)
    @assert length(is_eq) == length(cons) "The number of constraints must be the same as the number of equality conditions."
    vars = sorted_union(variables(objective), [variables(c) for c in cons]...)
    if !isempty(comm_gps)
        @assert issubset(union(comm_gps...), vars) "The commutative variables must be a subset of the variables."
    else
        push!(comm_gps, Set(vars))
    end
    @assert !(is_unipotent && is_projective) "The problem cannot be both unipotent and projective."
    return PolyOpt{T,obj_type}(objective, cons, is_eq, vars, comm_gps, is_unipotent, is_projective)
end

struct StatePolyOpt{T} <: OptimizationProblem
    objective::StatePolynomial{T}
    constraints::Vector{NCStatePolynomial{T}}
    is_equality::Vector{Bool}
    variables::Vector{Variable}
    comm_gps::Vector{Set{Variable}} # Set of variables that commutes with variables not in the set
    is_unipotent::Bool # square to 1. Examples: Pauli Operators, SWAP Operators
    is_projective::Bool # X^2 = X. Is projective.
end

function StatePolyOpt(objective::StatePolynomial{T}; constraints=Any[], is_equality=fill(false, length(constraints)), comm_gps=Vector{Variable}[], is_unipotent::Bool=false, is_projective::Bool=false) where {T}
    cons = collect(NCStatePolynomial{T}, constraints)
    is_eq = collect(Bool, is_equality)
    @assert length(is_eq) == length(cons) "The number of constraints must be the same as the number of equality conditions."
    vars = sorted_union(variables(objective), [variables(c) for c in cons]...)
    isempty(comm_gps) && push!(comm_gps, vars)
    @assert all([isempty(intersect(gp_a, gp_b)) for gp_a in comm_gps, gp_b in comm_gps if gp_a != gp_b]) "The commutative groups must be disjoint."
    @assert sorted_union(comm_gps...) == sort(vars) "The commutative variables groups must be equivalent to all the variables"
    @assert !(is_unipotent && is_projective) "The problem cannot be both unipotent and projective."
    return StatePolyOpt{T}(objective, cons, is_eq, vars, Set.(comm_gps), is_unipotent, is_projective)
end

function Base.show(io::IO, spolyopt::StatePolyOpt)
    cons_str = join(["$(string(con)) " * (iseq ? "= 0" : ">= 0") for (con, iseq) in zip(spolyopt.constraints, spolyopt.is_equality)], " \n")
    res_str = """
        obj: $(string(spolyopt.objective)) \n
        constraints: \n
            $(cons_str) \n
        variables:
            $(join(string.(spolyopt.variables)," ")) \n
        is_unipotent:
            $(spolyopt.is_unipotent) \n
        is_projective:
            $(spolyopt.is_projective) \n
    """
    print(io, res_str)
end
