using Test, NCTSSoS
using NCTSSoS.FastPolynomials: Variable, variables

@testset "PolyOpt Constructor" begin
    nvars = 10
    ncons = 3
    @ncpolyvar x[1:nvars]
    objective = 1.0 * sum(x .^ 2)
    constraints = [1.0 * sum(i .* x) for i in 1:ncons]

    @testset "Unconstrained" begin
        pop = PolyOpt(objective)

        @test pop.is_equality == Bool[]
        @test sort(pop.variables) == sort(x)
        @test pop.comm_gps == Set{Variable}[]
        @test !pop.is_unipotent
        @test !pop.is_projective

        comm_gps = [Set([x[1]])]

        pop = PolyOpt(objective; comm_gps=[Set([x[1]])], obj_type=TRACE)

        @test pop.comm_gps == [Set([x[1]])]
        @test pop isa PolyOpt{Float64,TRACE}
    end

    @testset "Constrainted Optimization Problem" begin
        pop = PolyOpt(objective; constraints=constraints)

        @test pop.constraints == constraints
        @test pop.is_equality == fill(false, ncons)

        pop = PolyOpt(objective; constraints=Set([constraints; sum(x)]))

        hash(pop.constraints[1])
        hash(sum(x))
         == hash(constraints)
        Set([constraints; sum(x)])

        constraints[1] == sum(x)

        @test length(pop.constraints) == ncons
        @test pop.is_equality == fill(false, ncons)

        is_equality = [isodd(i) ? true : false for i in 1:ncons]
        pop = PolyOpt(objective; constraints=constraints, is_equality=is_equality,is_unipotent=false,is_projective=true)
        @test pop.is_unipotent == false
        @test pop.is_projective == true
        @test pop.is_equality == is_equality
    end

    @testset "Invalid Input" begin
        @test_throws AssertionError PolyOpt(objective; is_equality= [true])
        @test_throws AssertionError PolyOpt(objective; constraints= constraints, is_equality=fill(true, ncons + 1), is_unipotent=false, is_projective=false)
        @test_throws AssertionError PolyOpt(objective; constraints= constraints, is_unipotent=true, is_projective=true)
        @test_throws AssertionError PolyOpt(x[1]*x[2]+x[3]*x[2])
        @ncpolyvar y[1:nvars]
        @test_throws AssertionError PolyOpt(objective; comm_gp=y)
    end
end


# using NCTSSoS: StatePolynomial, StateWord, coefficient
# @testset "StatePolyOpt Constructor" begin
#     @testset "Example 7.2.1" begin
#         @ncpolyvar x[1:2] y[1:2]
#         sp1 = sum([1.0, 1.0] .* map(a -> prod(ς.(a)), [[x[1] * y[2]], [y[2] * x[1]]]))
#         sp2 = sum([1.0, -1.0] .* map(a -> prod(ς.(a)), [[x[1] * y[1]], [x[2] * y[2]]]))
#         sp = sum([sp1 * sp1, sp2 * sp2])

#         sp1_sq = sum([1.0, 1.0, 1.0, 1.0] .* map(a -> prod(ς.(a)), [[x[1] * y[2], x[1] * y[2]], [y[2] * x[1], y[2] * x[1]], [x[1] * y[2], y[2] * x[1]], [y[2] * x[1], x[1] * y[2]]]))
#         sp2_sq = sum([1.0, -1.0, -1.0, 1.0] .* map(a -> prod(ς.(a)), [[x[1] * y[1], x[1] * y[1]], [x[1] * y[1], x[2] * y[2]], [x[2] * y[2], x[1] * y[1]], [x[2] * y[2], x[2] * y[2]]]))
#         true_obj = sum([sp1_sq, sp2_sq])

#         pop = StatePolyOpt(sp; is_unipotent=true, comm_gps=[x, y])
#         @test pop.objective ==  true_obj
#         @test pop.constraints == []
#         @test pop.is_equality == Bool[]
#         @test pop.is_unipotent == true 
#         @test pop.is_projective == false
#         @test pop.comm_gps == [Set(x),Set(y)]
#     end
#     @testset "Example 7.2.2" begin
#         @ncpolyvar x[1:3] y[1:3]
#         cov(a, b) = 1.0 * NCStateWord([x[a] * y[b]], one(x[1])) - 1.0 * (NCStateWord(monomial.([x[a]]), one(x[1])) * NCStateWord(monomial.([y[b]]), one(x[1])))
#         sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)

#         pop = StatePolyOpt(sp; is_unipotent=true,comm_gps= [x,y])
#         true_obj = sum([1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0] .* map(a -> prod(ς.(a)), ([[x[1] * y[1]], [x[1], y[1]], [x[1] * y[2]], [x[1], y[2]], [x[1] * y[3]], [x[1], y[3]], [x[2] * y[1]], [x[2], y[1]], [x[2] * y[2]], [x[2], y[2]], [x[2] * y[3]], [x[2], y[3]], [x[3] * y[1]], [x[3], y[1]], [x[3] * y[2]], [x[3], y[2]]])))
#         @test pop.objective == true_obj
#         @test pop.constraints == []
#         @test pop.is_unipotent == true
#         @test pop.is_equality == Bool[]
#         @test pop.comm_gps == [Set(x),Set(y)]
        
#         @ncpolyvar z[1:3]
#         @test_throws AssertionError StatePolyOpt(sp; comm_gps=[x, y, z])
#     end

#     @testset "Example 8.1.2" begin
#         @ncpolyvar A[1:3] B[1:3] 
#         J1 = sum(0.5 .* [coefficient(t) * ς(monomial(t)) for t in terms(1.0 * (A[1] + A[2] + A[3] + one(A[1]) * B[1] * (A[1] + A[2])))])
#         J2 = sum(0.5 .* [coefficient(t) * ς(monomial(t)) for t in terms(1.0 * (A[1] + A[2] - A[3] + one(A[1]) * B[2] * (A[1] - A[2])))]) + sum(0.5 .* [coefficient(t) * ς(monomial(t)) for t in terms(1.0 * (A[1] - A[2]) * (B[3] * A[1] - B[3] * A[2]))])
#         L = sum([4.0, 1.0, 1.0] .* [ς(monomial(v)) for v in [one(A[1]), A[1], A[2]]])

#         sp = sum([2.0 * J1 * J2, 2.0 * J1 * L, 2.0 * J2 * L, -1.0 * J1 * J1, -1.0 * J2 * J2, -1.0 * L * L])
#         pop = StatePolyOpt(sp; is_unipotent=true, comm_gps = [A,B])
#         @test pop.constraints == []
#         @test pop.is_unipotent == true
#         @test pop.is_equality == Bool[]
#         @test pop.comm_gps == [Set(A),Set(B)]
#     end
# end

