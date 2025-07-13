struct SOSProblem{T} <: OptimizationProblem
    model::GenericModel{T}
end

# Decompose the matrix into the form sum_j C_αj * g_j
# j: index of the constraint
# α: the monomial (JuMP variable)
function get_Cαj(basis_dict::Dict{GenericVariableRef{T},Int}, localizing_mtx::VectorConstraint{F,S,Shape}) where {T,F,S,Shape}
    dim = get_dim(localizing_mtx)
    cis = CartesianIndices((dim, dim))

    # basis idx, row, col
    dictionary_of_keys = Dict{Tuple{Int,Int,Int},T}()

    for (ci, cur_expr) in zip(cis, localizing_mtx.func)
        for (α, coeff) in cur_expr.terms
            dictionary_of_keys[(basis_dict[α], ci.I[1], ci.I[2])] = coeff
        end
    end

    return dictionary_of_keys
end

function get_Cαj_complex(basis_dict_real::Dict{GenericVariableRef{T},Int}, basis_dict_imag::Dict{GenericVariableRef{T},Int}, localizing_mtx::VectorConstraint{F,S,Shape}) where {T,F,S,Shape}
    dim = get_dim(localizing_mtx)
    cis = CartesianIndices((dim, dim))

    # basis idx, row, col
    dictionary_of_keys_real = Dict{Tuple{Int,Int,Int},Complex{T}}()
    dictionary_of_keys_imag = Dict{Tuple{Int,Int,Int},Complex{T}}()

    for (ci, cur_expr) in zip(cis, localizing_mtx.func)
        for (α, coeff) in cur_expr.terms
            if haskey(basis_dict_real, α)
                triple = (basis_dict_real[α], ci.I[1], ci.I[2])
                # if !haskey(dictionary_of_keys,triple)
                    dictionary_of_keys_real[triple] = coeff
                # else
                #     dictionary_of_keys[triple] += coeff
                # end
            elseif haskey(basis_dict_imag, α)
                triple = (basis_dict_imag[α], ci.I[1], ci.I[2])
                # if !haskey(dictionary_of_keys,triple)
                    dictionary_of_keys_imag[triple] = im*coeff
                # else
                    # dictionary_of_keys_imag[triple] += im*coeff
                # end
            else
                error("Unknown variable: $α")
            end
        end
    end
    return dictionary_of_keys_real, dictionary_of_keys_imag
end

"""
    sos_dualize(moment_problem::MomentProblem{T,M}) where {T,M} -> SOSProblem

Convert a moment problem into its dual SOS (Sum of Squares) problem formulation.

This function takes a moment problem and constructs the corresponding dual optimization
problem by introducing matrix variables for each constraint and formulating the dual
constraints that ensure polynomial equality.

# Arguments
- `moment_problem::MomentProblem{T,M}`: The primal moment problem to dualize

# Returns
- `SOSProblem`: The dual SOS problem with matrix variables and constraints

# Details
The dualization process involves:
1. Creating matrix variables (G_j) for each constraint, either in symmetric matrix space
   or positive semidefinite cone depending on the constraint type
2. Introducing a scalar variable `b` to bound the minimum value of the primal problem
3. Setting up polynomial equality constraints by matching coefficients of monomials
4. Handling symmetrization of the monomial basis to ensure proper polynomial comparison

The resulting dual problem maximizes `b` subject to the constraint that the sum of
matrix variables weighted by coefficient matrices equals the objective polynomial.
"""
function sos_dualize(moment_problem::MomentProblem{T,M}) where {T,M}
    dual_model = GenericModel{T}()

    # Initialize Gj as variables
    dual_variables = map(constraint_object.(moment_problem.constraints)) do cons
        G_dim = get_dim(cons)
        @variable(dual_model, [1:G_dim, 1:G_dim] in ((cons.set isa MOI.Zeros) ? SymmetricMatrixSpace() : (moment_problem.complex_coeff ? HermitianPSDCone() : PSDCone())))
    end

    if moment_problem.complex_coeff
        # b: to bound the minimum value of the primal problem
        @variable(dual_model, b)
        @objective(dual_model, Max, b)

        primal_objective_terms = objective_function(moment_problem.model).terms

        # NOTE: objective is Symmetric, hence when comparing polynomials, we need to canonicalize them first
        unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))

        symmetric_basis = unsymmetrized_basis 

        # JuMP variables corresponding to symmetric_basis
        symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)

        # specify constraints
        fα_constraints_real = [GenericAffExpr{Complex{T},VariableRef}(get(primal_objective_terms, α.terms.keys[1], zero(Complex{T}))) for α in symmetric_variables]

        # fα_constraints_imag = [GenericAffExpr{Complex{T},VariableRef}(get(primal_objective_terms, α.terms.keys[2], zero(Complex{T}))) for α in symmetric_variables]

        symmetrized_α2cons_dict = Dict(zip(unsymmetrized_basis, map(x -> searchsortedfirst(symmetric_basis, canonicalize(x, moment_problem.sa)), unsymmetrized_basis)))

        unsymmetrized_basis_vals_dict_real = Dict(zip(getindex.(getfield.(getfield.(getindex.(Ref(moment_problem.monomap), unsymmetrized_basis), Ref(:terms)), Ref(:keys)), Ref(1)), 1:length(unsymmetrized_basis)))

        unsymmetrized_basis_vals_dict_imag = Dict(zip(getindex.(getfield.(getfield.(getindex.(Ref(moment_problem.monomap), unsymmetrized_basis), Ref(:terms)), Ref(:keys)), Ref(2)), 1:length(unsymmetrized_basis)))

        add_to_expression!(fα_constraints_real[1], -one(Complex{T}), b)

        for (i, sdp_constraint) in enumerate(moment_problem.constraints)
            Cαjs_real, Cαjs_imag = get_Cαj_complex(unsymmetrized_basis_vals_dict_real, unsymmetrized_basis_vals_dict_imag, constraint_object(sdp_constraint))
            @show Cαjs_real, Cαjs_imag
            for (ky, coef) in Cαjs_real
                @show ky coef, real(dual_variables[i][ky[2], ky[3]])
                add_to_expression!(fα_constraints_real[symmetrized_α2cons_dict[unsymmetrized_basis[ky[1]]]], -coef, real(dual_variables[i][ky[2], ky[3]]))
            end

            for (ky, coef) in Cαjs_imag
                @show ky coef, imag(dual_variables[i][ky[2], ky[3]])
                add_to_expression!(fα_constraints_real[symmetrized_α2cons_dict[unsymmetrized_basis[ky[1]]]], -coef, imag(dual_variables[i][ky[2], ky[3]]))
            end
        end
        @constraint(dual_model, fα_constraints_real .== 0)
        # @constraint(dual_model, fα_constraints_imag .== 0)

        return SOSProblem(dual_model)
    else
        # b: to bound the minimum value of the primal problem
        @variable(dual_model, b)
        @objective(dual_model, Max, b)

        primal_objective_terms = objective_function(moment_problem.model).terms

        # NOTE: objective is Symmetric, hence when comparing polynomials, we need to canonicalize them first
        unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))

        symmetric_basis = sorted_unique(canonicalize.(unsymmetrized_basis, Ref(moment_problem.sa)))

        # JuMP variables corresponding to symmetric_basis
        symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)


        fα_constraints = [GenericAffExpr{T,VariableRef}(get(primal_objective_terms, α, zero(T))) for α in symmetric_variables]

        symmetrized_α2cons_dict = Dict(zip(unsymmetrized_basis, map(x -> searchsortedfirst(symmetric_basis, canonicalize(x, moment_problem.sa)), unsymmetrized_basis)))

        unsymmetrized_basis_vals_dict = Dict(zip(getindex.(Ref(moment_problem.monomap), unsymmetrized_basis), 1:length(unsymmetrized_basis)))

        add_to_expression!(fα_constraints[1], -one(T), b)

        for (i, sdp_constraint) in enumerate(moment_problem.constraints)
            Cαjs = get_Cαj(unsymmetrized_basis_vals_dict, constraint_object(sdp_constraint))
            for (ky, coef) in Cαjs
                add_to_expression!(fα_constraints[symmetrized_α2cons_dict[unsymmetrized_basis[ky[1]]]], -coef, dual_variables[i][ky[2], ky[3]])
            end
        end
        @constraint(dual_model, fα_constraints .== 0)

        return SOSProblem(dual_model)
    end
end
