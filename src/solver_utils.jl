function get_dim(cons::VectorConstraint)
    if typeof(JuMP.shape(cons)) in (JuMP.HermitianMatrixShape, JuMP.SquareMatrixShape)
        return JuMP.shape(cons).side_dimension
    elseif typeof(JuMP.shape(cons)) <: JuMP.ArrayShape
        return JuMP.shape(cons).dims[1]
    else
        return error("Unsupported constraint shape $(typeof(JuMP.shape(cons)))")
    end
end