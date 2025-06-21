# T: type of the coefficients
# monomap: map from monomials in DynamicPolynomials to variables in JuMP
struct MomentProblem{T,M,CR<:ConstraintRef} <: OptimizationProblem
    model::GenericModel{T}
    constraints::Vector{CR}
    monomap::Dict{M,GenericVariableRef{T}}  # TODO: maybe refactor.
    sa::SimplifyAlgorithm
end

function substitute_variables(poly::P, monomap::Dict{M,GenericVariableRef{T}}) where {T,P<:AbstractPolynomial{T},M}
    mapreduce(+, zip(coefficients(poly), monomials(poly))) do (coef, mono)
        coef * monomap[expval(mono)]
    end
end

function get_mom_matrix(mom_problem::MomentProblem)
    _, mom_loc = findmax(get_dim, constraint_object.(mom_problem.constraints))
    return value.(mom_problem.constraints[mom_loc])
end

# cliques_cons: groups constraints according to cliques,
# global_cons: constraints that are not in any single clique
# cliques_term_sparsities: each clique, each obj/constraint, each ts_clique, each basis needed to index moment matrix
function moment_relax(pop::PolyOpt{P}, corr_sparsity::CorrelativeSparsity, cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}) where {T,P<:AbstractPolynomial{T},M}
    # NOTE: objective and constraints may have integer coefficients, but popular JuMP solvers does not support integer coefficients
    # left type here to support BigFloat model for higher precision
    model = GenericModel{T}()

    sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, is_unipotent=pop.is_unipotent, is_projective=pop.is_projective)
    # the union of clique_total_basis
    total_basis = sorted_union(map(zip(corr_sparsity.clq_cons, cliques_term_sparsities)) do (cons_idx, term_sparsities)
        union(vec(reduce(vcat, [
            map(monomials(poly)) do m
                expval(simplify(neat_dot(rol_idx, m * col_idx), sa))
            end
            for (poly, term_sparsity) in zip([one(pop.objective); corr_sparsity.cons[cons_idx]], term_sparsities) for basis in term_sparsity.block_bases for rol_idx in basis for col_idx in basis
        ])))
    end...)

    for (cons_idx, term_sparsities) in zip(corr_sparsity.clq_cons, cliques_term_sparsities)
        for (poly, term_sparsity) in zip([one(pop.objective); corr_sparsity.cons[cons_idx]], term_sparsities)
            @show poly
            for basis in term_sparsity.block_bases
                @show basis 
                for rol_idx in eachindex(basis) 
                    for col_idx in eachindex(basis) 
                        map(monomials(poly)) do m
                            @show rol_idx, col_idx, expval(simplify(neat_dot(basis[rol_idx], m * basis[col_idx]), sa))
                        end
                    end
                end
            end
        end
    end

    # map the monomials to JuMP variables, the first variable must be 1
    @variable(model, y[1:length(total_basis)])
    @constraint(model, first(y) == 1)
    monomap = Dict(zip(total_basis, y))

    @show keys(monomap)
    for tss in cliques_term_sparsities
        @show "new clique"
        for ts in tss
            @show "new poly"
            for bb in ts.block_bases
                @show bb
            end
        end
    end

    constraint_matrices =
        [mapreduce(vcat, zip(cliques_term_sparsities, corr_sparsity.clq_cons)) do (term_sparsities, cons_idx)
                mapreduce(vcat, zip(term_sparsities, [one(pop.objective), corr_sparsity.cons[cons_idx]...])) do (term_sparsity, poly)
                    map(term_sparsity.block_bases) do ts_sub_basis
                        constrain_moment_matrix!(
                            model,
                            poly,
                            ts_sub_basis,
                            monomap,
                            poly in pop.eq_constraints ? Zeros() : PSDCone(), sa)
                    end
                end
            end
            map(corr_sparsity.global_cons) do global_con
                constrain_moment_matrix!(
                    model,
                    corr_sparsity.cons[global_con],
                    [one(pop.objective)],
                    monomap,
                    global_con <= length(pop.eq_constraints) ? Zeros() : PSDCone(),
                    sa
                )
            end]

    @objective(model, Min, substitute_variables(mapreduce(p -> p[1] * simplify(p[2], sa), +, terms(symmetric_canonicalize(pop.objective, sa)); init=zero(pop.objective)), monomap))

    return MomentProblem(model, constraint_matrices, monomap, sa)
end

function constrain_moment_matrix!(
    model::GenericModel{T},
    poly::P,
    local_basis::Vector{M1}, # M2 should be expval(M1)
    monomap::Dict{M2,GenericVariableRef{T}},
    cone, # FIXME: which type should I use?
    sa::SimplifyAlgorithm
) where {T,P<:AbstractPolynomial{T},M1,M2}
    moment_mtx = [
        substitute_variables(sum([coef * simplify(neat_dot(row_idx, mono * col_idx), sa) for (coef, mono) in zip(coefficients(poly), monomials(poly))]), monomap) for
        row_idx in local_basis, col_idx in local_basis
    ]
    return @constraint(model, moment_mtx in cone)
end
