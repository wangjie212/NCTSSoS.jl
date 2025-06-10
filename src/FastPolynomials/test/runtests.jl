using NCTSSoS.FastPolynomials
using Test

@testset "FastPolynomials.jl" begin
    include("arithmetic.jl")
    include("compare.jl")
    include("monomials.jl")
    include("polynomial.jl")
    include("statepolynomial.jl")
    include("utils.jl")
    include("variables.jl")
end
