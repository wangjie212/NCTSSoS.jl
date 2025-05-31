using Test, NCTSSoS
using JuMP, DynamicPolynomials
using NCTSSoS: get_state_basis, neat_dot, NCStateWord, NCStatePolynomial, constrain_moment_matrix!, expval, substitute_variables, NCStateTerm, moment_relax
using Clarabel, COSMO
using NCTSSoS: sos_dualize
using NCTSSoS: correlative_sparsity, iterate_term_sparse_supp, sorted_union, MinimalChordal, NoElimination

@testset "State Polynomial Opt 7.2.0" begin
    @ncpolyvar x[1:2] y[1:2]
    sp = NCStatePolynomial(map(a -> a[1]*NCStateWord([a[2]],a[3]), zip([-1.0, -1.0, -1.0, 1.0], [x[1] * y[1], x[1] * y[2], x[2] * y[1], x[2] * y[2]],fill(one(x[1]),4))))
    spop = StatePolyOpt(sp; is_unipotent=true, comm_gps=[x, y])

    d = 1
    cr = correlative_sparsity(spop, d, NoElimination())

    cliques_objective = [reduce(+, [issubset(effective_variables(t.ncstate_word), clique) ? t : zero(t) for t in terms(spop.objective)]) for clique in cr.cliques]

    initial_activated_supp = [sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, spop.constraints[cons_idx]; init=typeof(monomials(spop.objective)[1])[]))
                              for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, cr.cliques_cons, cr.cliques_idcs_bases)]

    cliques_term_sparsities = map(zip(initial_activated_supp, cr.cliques_cons, cr.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
        [iterate_term_sparse_supp(activated_supp, poly, basis, NoElimination()) for (poly, basis) in zip([one(spop.objective); spop.constraints[cons_idx]], idcs_bases)]
    end

    mom_problem = moment_relax(spop, cr.cliques_cons, cr.global_cons, cliques_term_sparsities)

    # FIXME: cannot find some key
    sos_problem = sos_dualize(mom_problem)

    set_optimizer(mom_problem.model, Clarabel.Optimizer)
    optimize!(mom_problem.model)
    @test isapprox(objective_value(mom_problem.model), -2.8284271321623202, atol=1e-5)
    @test is_solved_and_feasible(mom_problem.model)

    set_optimizer(sos_problem.model, Clarabel.Optimizer)
    optimize!(sos_problem.model)
    @test isapprox(objective_value(sos_problem.model), -2.8284271321623202, atol=1e-5)
    @test is_solved_and_feasible(sos_problem.model)
end

@testset "State Polynomial Opt 7.2.1" begin
    @ncpolyvar x[1:2] y[1:2]
    sp1 = NCStatePolynomial([coef*sw for (coef,sw) in zip([1.0,1.0],NCStateWord.([[x[1] * y[2]], [x[2] * y[1]]], Ref(one(x[1]))))])
    sp2 = NCStatePolynomial([coef*sw for (coef,sw) in zip([1.0, -1.0],NCStateWord.([[x[1] * y[1]], [x[2] * y[2]]], Ref(one(x[1]))))])
    words = [one(x[1]), one(x[1])]
    sp = -1.0 * sp1 * sp1  + (-1.0 * sp2 * sp2)

    spop = StatePolyOpt(sp; is_unipotent=true, comm_gps=[x, y])

    d = 3
    cr = correlative_sparsity(spop, d, NoElimination())

    cliques_objective = [reduce(+, [issubset(effective_variables(t.ncstate_word), clique) ? t : zero(t) for t in terms(spop.objective)]) for clique in cr.cliques]

    initial_activated_supp = [sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, spop.constraints[cons_idx]; init=typeof(monomials(spop.objective)[1])[]))
                              for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, cr.cliques_cons, cr.cliques_idcs_bases)]

    cliques_term_sparsities = map(zip(initial_activated_supp, cr.cliques_cons, cr.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
        [iterate_term_sparse_supp(activated_supp, poly, basis, NoElimination()) for (poly, basis) in zip([one(spop.objective); spop.constraints[cons_idx]], idcs_bases)]
    end

    mom_problem = moment_relax(spop, cr.cliques_cons, cr.global_cons, cliques_term_sparsities)

    set_optimizer(mom_problem.model, COSMO.Optimizer)
    optimize!(mom_problem.model)
    @test isapprox(objective_value(mom_problem.model), -4.0, atol=1e-5)
    @test is_solved_and_feasible(mom_problem.model)

    sos_problem = sos_dualize(mom_problem)
    set_optimizer(sos_problem.model, COSMO.Optimizer)
    optimize!(sos_problem.model)
    @test isapprox(objective_value(sos_problem.model), -4.0, atol=1e-5)
    @test is_solved_and_feasible(sos_problem.model)
end

@testset "State Polynomial Opt 7.2.2" begin
    @ncpolyvar x[1:6]
    sp = NCStatePolynomial(map(a->a[1]*NCStateWord(monomial.(a[2]),one(x[1])),zip(Float64.([-1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1]),[[x[1] * x[4]], [x[1], x[4]], [x[1] * x[5]], [x[1], x[5]], [x[1] * x[6]], [x[1], x[6]], [x[2] * x[4]], [x[2], x[4]], [x[2] * x[5]], [x[2], x[5]], [x[2] * x[6]], [x[2], x[6]], [x[3] * x[4]], [x[3], x[4]], [x[3] * x[5]], [x[3], x[5]]])))

    spop = StatePolyOpt(sp; is_unipotent=true, comm_gps=[x[1:3], x[4:6]])

    d = 2
    cr = correlative_sparsity(spop, d, NoElimination())

    cliques_objective = [reduce(+, [issubset(effective_variables(t.ncstate_word), clique) ? t : zero(t) for t in terms(spop.objective)]) for clique in cr.cliques]

    initial_activated_supp = [sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, spop.constraints[cons_idx]; init=typeof(monomials(spop.objective)[1])[]))
                              for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, cr.cliques_cons, cr.cliques_idcs_bases)]

    cliques_term_sparsities = map(zip(initial_activated_supp, cr.cliques_cons, cr.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
        [iterate_term_sparse_supp(activated_supp, poly, basis, NoElimination()) for (poly, basis) in zip([one(spop.objective); spop.constraints[cons_idx]], idcs_bases)]
    end

    mom_problem = moment_relax(spop, cr.cliques_cons, cr.global_cons, cliques_term_sparsities)
    
    sos_problem = sos_dualize(mom_problem)
    set_optimizer(sos_problem.model, optimizer_with_attributes(COSMO.Optimizer,"eps_rel"=> 1e-8))
    optimize!(sos_problem.model)
    @test isapprox(objective_value(sos_problem.model), -5.0, atol=1e-3)
end


@testset "Constrain Moment matrix" begin
    @ncpolyvar x[1:2]

    basis = get_state_basis(x,1,identity)

    sp = NCStatePolynomial(map(a -> a[1]* a[2], zip([1.0, 2.0, 3.0], NCStateWord.(map(x -> monomial.(x), [[x[1], x[2]], [x[1]], [x[2]]]),Ref(one(x[1]))))))
    nc_words = monomial.([one(x[1]), x[1], x[2]])
    ncsp = NCStatePolynomial(map(a -> a[1] * a[2], zip([1.0, 2.0, 3.0], NCStateWord.(map(x -> monomial.(x), [[x[1], x[2]], [x[1]], [x[2]]]),nc_words))))
    poly = one(ncsp)

    total_basis = sort(unique([expval(neat_dot(a,b)) for a in basis for b in basis]))

    model = GenericModel{Float64}()
    @variable(model, y[1:length(total_basis)])
    wordmap = Dict(zip(total_basis,y))

    ncterms = map(a -> a[1]* NCStateWord(monomial.(a[2]), a[3]), zip([1.0, 2.0, 3.0], [[x[1], x[2]], [x[1]], [x[2]]], nc_words))
    @test sort(terms(ncsp)) == sort(ncterms)

    @test substitute_variables(ncsp, wordmap) == 1.0 * y[4] + 3.0 * y[3] + 2.0 * y[6]

    true_mom_mtx = expval.([neat_dot(a,b) for a in basis, b in basis])
    mom_mtx_cons = constrain_moment_matrix!(model, one(ncsp), basis, wordmap, PSDCone(), identity)
    mom_mtx = constraint_object(mom_mtx_cons)
    @test reshape(mom_mtx.func, 5, 5) == AffExpr[y[1] y[2] y[5] y[2] y[5]; y[2] y[3] y[4] y[3] y[4]; y[5] y[4] y[6] y[4] y[6]; y[2] y[3] y[4] y[8] y[7]; y[5] y[4] y[6] y[9] y[10]]
end


