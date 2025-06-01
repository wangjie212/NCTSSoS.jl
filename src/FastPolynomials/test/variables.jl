using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: Variable, COMPLEX, REAL
using NCTSSoS.FastPolynomials: polyarrayvar, buildpolyvar, buildpolyvars

@testset "Variable" begin
    @testset "Creation by Macros" begin
        my_vars = @ncpolyvar x y z
        @test my_vars isa Tuple
        @test x isa Variable
        @test y isa Variable
        @test x.name == "x"
        @test y.name == "y"

        @ncpolyvar x[1:10]
        @test x isa Vector{Variable}
        @test length(x) == 10
        @test x[1].name == "x[1]"

        @ncpolyvar x[1:10, 1:5]
        @test x isa Matrix{Variable}
        @test size(x) == (10, 5)
    end

    @testset "ComplexKind Conversion" begin
        @ncpolyvar x_real
        x_complex = Variable(x_real, COMPLEX)
        @test x_complex isa Variable
        @test x_complex.name == "x_real"
        @test x_complex.kind == COMPLEX
    end

    @testset "Utils: polyarrayvar" begin
        x = polyarrayvar(REAL, "x", 1:3, 1:2)
        @test x isa Matrix{Variable}
        @test size(x) == (3, 2)
    end

    @testset "Utils: buildpolyvar" begin
        @test_throws ErrorException buildpolyvar(Expr(:call, :x, +), REAL)
    end
end
