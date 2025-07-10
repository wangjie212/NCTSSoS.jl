using NCTSSoS, Test

@testset "NCTSSoS.jl" begin
    include("fastpoly_test/runtests.jl")
    include("pop.jl")
    include("sparse.jl")
    include("solver_utils.jl")
    include("moment_solver.jl")
    include("sos_solver.jl")
    include("interface.jl")
    include("state_poly_opt.jl")
    include("Aqua.jl")
    include("Doctest.jl")
    include("ExplicitImports.jl")
    include("past_bugs.jl")
end
