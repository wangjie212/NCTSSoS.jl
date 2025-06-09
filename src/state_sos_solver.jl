function sos_dualize(moment_problem::StateMomentProblem{T}) where {T}
    dual_model = GenericModel{T}()

    # Initialize Gj as variables
    dual_variables = map(constraint_object.(moment_problem.constraints)) do cons
        # FIXME: HEDIOUS
        G_dim = get_dim(cons)
        @variable(dual_model, [1:G_dim, 1:G_dim] in ((cons.set isa MOI.Zeros) ? SymmetricMatrixSpace() : PSDCone()))
    end

    # b: to bound the minimum value of the primal problem
    @variable(dual_model, b)
    @objective(dual_model, Max, b)

    primal_objective_terms = objective_function(moment_problem.model).terms

    unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))

    symmetric_basis = unsymmetrized_basis

    # JuMP variables corresponding to symmetric_basis
    symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)

    # specify constraints
    fα_constraints = [AffExpr(get(primal_objective_terms, α, zero(T))) for α in symmetric_variables]

    symmetrized_α2cons_dict = Dict(zip(unsymmetrized_basis, map(x -> searchsortedfirst(symmetric_basis, prod(moment_problem.reduce_func(symmetric_canonicalize(x)))), unsymmetrized_basis)))

    unsymmetrized_basis_vals = getindex.(Ref(moment_problem.monomap), unsymmetrized_basis)

    add_to_expression!(fα_constraints[1], -1.0, b)

    for (i, sdp_constraint) in enumerate(moment_problem.constraints)
        Cαj = get_Cαj(unsymmetrized_basis_vals, constraint_object(sdp_constraint))
        for (k, α) in enumerate(unsymmetrized_basis)
            for jj in 1:size(Cαj[k], 2)
                for ii in nzrange(Cαj[k], jj)
                    add_to_expression!(fα_constraints[symmetrized_α2cons_dict[α]], -Cαj[k].nzval[ii], dual_variables[i][jj, Cαj[k].rowval[ii]])
                end
            end
        end
    end
    @constraint(dual_model, fα_constraints .== 0)

    return SOSProblem(dual_model)
end
