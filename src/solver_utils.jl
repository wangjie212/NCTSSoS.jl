function get_dim(cons::VectorConstraint)
    cons_shape_type = typeof(JuMP.shape(cons))
    if cons_shape_type == JuMP.HermitianMatrixShape
        return JuMP.shape(cons).side_dimension
    elseif cons_shape_type == JuMP.SquareMatrixShape
        return JuMP.shape(cons).side_dimension
    elseif cons_shape_type <: JuMP.ArrayShape
        return JuMP.shape(cons).dims[1]
    else
        return error("Unsupported constraint shape $(typeof(JuMP.shape(cons)))")
    end
end

get_coef_type(cons::VectorConstraint) = eltype(cons.func).parameters[1]
