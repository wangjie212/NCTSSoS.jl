using Test, NCTSSoS
using NCTSSoS: StatePolyOpt, NCStateWord
using NCTSSoS: AsIsElimination
using NCTSSoS: sorted_union, symmetric_canonicalize, neat_dot, iterate_term_sparse_supp
using NCTSSoS: get_correlative_graph, correlative_sparsity, moment_relax, sos_dualize
using DynamicPolynomials
using Graphs
using COSMO
using JuMP
using Clarabel
using MosekTools

@testset "Correlative Sparsity CHSH" begin
    @ncpolyvar x[1:2] y[1:2]
    sp = -1.0 *  ς(x[1]*y[1]) * one(x[1]) - 1.0 * ς(x[1]*y[2]) * one(x[1]) - (1.0 * ς(x[2]*y[1]) * one(x[1]) ) + 1.0 * ς(x[2]*y[2]) * one(x[1])
    spop = StatePolyOpt(sp; is_unipotent=true, comm_gps=[x, y])

    d = 3
    cg = get_correlative_graph(spop.variables, spop.objective, spop.constraints, d)
    @test cg.fadjlist == [[3,4],[3,4],[1,2],[1,2]]

    cr = correlative_sparsity(spop, d, NoElimination())
    @test cr.cliques == [[x;y]]
    @test cr.cliques_cons == [Int64[]]
    @test cr.global_cons == Int64[]

    cliques_objective = [reduce(+, [issubset(effective_variables(t.ncstate_word), clique) ? t : zero(t) for t in terms(spop.objective)]) for clique in cr.cliques]

    initial_activated_supp = [sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, spop.constraints[cons_idx]; init=typeof(monomials(spop.objective)[1])[]))
                              for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, cr.cliques_cons, cr.cliques_idcs_bases)]

    cliques_term_sparsities = map(zip(initial_activated_supp, cr.cliques_cons, cr.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
        [iterate_term_sparse_supp(activated_supp, poly, basis, AsIsElimination()) for (poly, basis) in zip([one(spop.objective); spop.constraints[cons_idx]], idcs_bases)]
    end

    mom_problem = moment_relax(spop, cr.cliques_cons, cr.global_cons, cliques_term_sparsities)
    set_optimizer(mom_problem.model, Clarabel.Optimizer)
    optimize!(mom_problem.model)
    @test isapprox(objective_value(mom_problem.model), -2.828, atol=1e-3)
    @test is_solved_and_feasible(mom_problem.model)

    # why am I missing a term?
    # sos_problem = sos_dualize(mom_problem)
    # set_optimizer(sos_problem.model, Clarabel.Optimizer)
    # optimize!(sos_problem.model)
    # @test is_solved_and_feasible(sos_problem.model)
    # @test isapprox(objective_value(sos_problem.model), -2.828, atol=1e-3)
end

@testset "Correlative Sparsity" begin
    @ncpolyvar x[1:2] y[1:2]
    sp1 = NCStatePolynomial([coef*sw for (coef,sw) in zip([1.0,1.0],NCStateWord.([[x[1] * y[2]], [x[2] * y[1]]], Ref(one(x[1]))))])
    sp2 = NCStatePolynomial([coef*sw for (coef,sw) in zip([1.0, -1.0],NCStateWord.([[x[1] * y[1]], [x[2] * y[2]]],Ref(one(x[1]))))])
    words = [one(x[1]), one(x[1])]
    sp = -1.0 * sp1 * sp1  + (-1.0 * sp2 * sp2)

    spop = StatePolyOpt(sp; is_unipotent=true, comm_gps=[x, y])

    using DynamicPolynomials: degree
    using NCTSSoS: symmetric_canonicalize, monomials
    cur_reducer = reducer(spop)
    nc_basis = filter(a -> degree(a) == 3, unique!(map(a -> prod(symmetric_canonicalize.(a)), (cur_reducer.(monomials([x; y], 3))))))
    nc_basis = filter(a -> degree(a) == 3, unique!(cur_reducer.(monomials([x; y], 3))))

    d = 3

    cg = get_correlative_graph(spop.variables, spop.objective, spop.constraints, d)
    @test cg.fadjlist == [[2,3,4],[1,3,4],[1,2,4],[1,2,3]]


    cr = correlative_sparsity(spop, d, NoElimination())
    @test cr.cliques == [[x;y]]
    @test cr.cliques_cons == [Int64[]]
    @test cr.global_cons == Int64[]

    cliques_objective = [reduce(+, [issubset(effective_variables(t.ncstate_word), clique) ? t : zero(t) for t in terms(spop.objective)]) for clique in cr.cliques]

    initial_activated_supp = [sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, spop.constraints[cons_idx]; init=typeof(monomials(spop.objective)[1])[]))
                              for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, cr.cliques_cons, cr.cliques_idcs_bases)]

    # state word basis is incorrect! should use trace basis but why?
    cliques_term_sparsities = map(zip(initial_activated_supp, cr.cliques_cons, cr.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
        [iterate_term_sparse_supp(activated_supp, poly, basis, AsIsElimination()) for (poly, basis) in zip([one(spop.objective); spop.constraints[cons_idx]], idcs_bases)]
    end
    # which variable group according to correlative sparsity
    # is it a moment matrix or localizing matrix?
    for el in cliques_term_sparsities[1][1].block_bases
        @show el[1]
    end

    mom_problem = moment_relax(spop, cr.cliques_cons, cr.global_cons, cliques_term_sparsities)
    set_optimizer(mom_problem.model, COSMO.Optimizer)
    optimize!(mom_problem.model)
    objective_value(mom_problem.model)
    @test isapprox(objective_value(mom_problem.model), -4.0, atol=1e-4)
    @test is_solved_and_feasible(mom_problem.model)

    # here too
    # sos_problem = sos_dualize(mom_problem)
    # set_optimizer(sos_problem.model, COSMO.Optimizer)
    # optimize!(sos_problem.model)
    # @test isapprox(objective_value(sos_problem.model), -4.0, atol=1e-4)
    # @test is_solved_and_feasible(sos_problem.model)
end