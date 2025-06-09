get_dim(cons::VectorConstraint) = cons.set isa MOI.PositiveSemidefiniteConeSquare ? JuMP.shape(cons).side_dimension : JuMP.shape(cons).dims[1]

function reducer(pop::PolyOpt)
    function (x)
        cxs = _comm(x, pop.comm_gps)
        return pop.is_unipotent ? _unipotent.(cxs) : (pop.is_projective ? _projective.(cxs) : cxs)
    end
end

function reducer(spop::StatePolyOpt)
    function (x)
        cxs = _comm(x, spop.comm_gps)
        return spop.is_unipotent ? _unipotent.(cxs) : (spop.is_projective ? _projective.(cxs) : cxs)
    end
end

function FastPolynomials._comm(ncsw::NCStateWord, comm_gps::Vector{Set{Variable}})
    [NCStateWord(prod.(_comm.(ncsw.sw, Ref(comm_gps))), prod(_comm(ncsw.nc_word, comm_gps)))]
end
