using FastPolynomials
using Test

@testset "FastPolynomials.jl" begin
    # Write your tests here.
    include("variables.jl")
    include("monomials.jl")
    include("arithmetic.jl")
    include("compare.jl")
    include("polynomial.jl")
end
