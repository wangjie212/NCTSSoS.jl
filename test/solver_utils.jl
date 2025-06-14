using Test, NCTSSoS
using JuMP
using NCTSSoS.FastPolynomials
using NCTSSoS: get_dim, reducer, sorted_unique

using NCTSSoS.FastPolynomials: Monomial, Polynomial, sorted_union, get_basis

@testset "Utilities" begin
    @ncpolyvar x y z
    @testset "VectorConstraint Diexport @ncpolyvar
m" begin
        model = Model()
        n = 5
        var1 = @variable(model, [1:n, 1:n])
        var2 = @variable(model, [1:2*n, 1:2*n])

        cons1 = @constraint(model, var1 in PSDCone())
        cons2 = @constraint(model, var2 in Zeros())

        @test get_dim(constraint_object(cons1)) == n
        @test get_dim(constraint_object(cons2)) == 2 * n
    end

    @testset "reducer" begin
        obj =1.0 * x * y + 2.0 * y * z

        basis = get_basis([x,y,z],3)

        pop= PolyOpt(obj; is_unipotent=true)
        reducer_func = reducer(pop)
        @test prod(reducer_func(y*x^2*y)) == one(Monomial) 

        pop = PolyOpt(obj; is_projective=true)
        reducer_func = reducer(pop)
        @test prod(reducer_func(y*x^2*y)) == y*x*y

        pop = PolyOpt(obj; comm_gps = [[x], [y,z]], is_unipotent=true)
        reducer_func = reducer(pop)
        @test prod(reducer_func(y*x^2*y)) == one(y)

        @test sorted_unique(map(prod, reducer_func.(basis))) == 
        sort([one(x * y), Monomial([z], [1]), Monomial([y], [1]), Monomial([x], [1]), z * y, x * y, y * z, x * z, x * z * y, z * y * z, y * z * y, x * y * z])

        pop = PolyOpt(obj; comm_gps=[[x], [y, z]], is_projective=true)
        reducer_func = reducer(pop)
        @test prod(reducer_func(y*x^2*y)) == x*y

        @test sorted_unique(map(prod, reducer_func.(basis))) == sort([one(x * y * z), Monomial([z], [1]), Monomial([y], [1]), Monomial([x], [1]), z * y, x * z, x * y, y * z, x * z * y, z * y * z, x * y * z, y * z * y])
    end
end
