get_dim(cons::VectorConstraint) = cons.set isa MOI.PositiveSemidefiniteConeSquare ? JuMP.shape(cons).side_dimension : JuMP.shape(cons).dims[1]

get_coef_type(cons::VectorConstraint) = eltype(cons.func).parameters[1]
