abstract type OptimizationProblem end

@enum Objective EIGEN TRACE
# T: type of the coefficients, currently removed to reduce redundence
# P: type of the polynomial, either `Polynomial{T}` or `NCStatePolynomial{T}`
# OBJ: type of the objective function, one of `Objective`
struct PolyOpt{P,OBJ} <: OptimizationProblem
    objective::P
    eq_constraints::Vector{P} # NOTE: assuming constraints are all simplified using comm_gp, is_unipotent, and is_projective
    ineq_constraints::Vector{P}
    variables::Vector{Variable}
    comm_gps::Vector{Vector{Variable}} # Vectors of Set of variables that commutes with variables not in the same set
    is_unipotent::Bool # square to 1. Examples: Pauli Operators, SWAP Operators
    is_projective::Bool # X^2 = X. Is projective.
end

function PolyOpt(objective::P; eq_constraints=Any[], ineq_constraints=Any[], comm_gps=Vector{Variable}[], is_unipotent::Bool=false, is_projective::Bool=false, obj_type::Objective=EIGEN) where {T,P<:AbstractPolynomial{T}}
    @assert !(T <: Integer) "The polynomial coefficients can not be integers (not supported by JuMP solvers)."
    eq_cons = unique!(collect(P, eq_constraints))
    ineq_cons = unique!(collect(P, ineq_constraints))
    vars = sorted_union(variables(objective), variables.(eq_cons)..., variables.(ineq_cons)...)
    if !isempty(comm_gps)
        @assert all([isempty(intersect(gp_a, gp_b)) for gp_a in comm_gps, gp_b in comm_gps if gp_a != gp_b]) "The commutative groups must be disjoint."
        @assert issubset(union(comm_gps...), vars) "The commutative variables must be a subset of the variables."
    else
        push!(comm_gps, vars)
    end
    @assert !(is_unipotent && is_projective) "The problem cannot be both unipotent and projective."
    return PolyOpt{P,obj_type}(objective, eq_cons, ineq_cons, vars, comm_gps, is_unipotent, is_projective)
end

function Base.show(io::IO, pop::PolyOpt)
    cons_str(cons::Vector{P}, iseq::Bool) where {P} = join(["$(string(c)) " * (iseq ? "= 0" : ">= 0") for c in cons], " \n")
    res_str = """
        obj: \n
            $(string(pop.objective)) \n
        constraints: \n
            $(cons_str(pop.eq_constraints,true)) \n
            $(cons_str(pop.ineq_constraints,false)) \n
        variables:
            $(join(string.(pop.variables)," ")) \n
        is_unipotent:
            $(pop.is_unipotent) \n
        is_projective:
            $(pop.is_projective) \n
    """
    print(io, res_str)
end
