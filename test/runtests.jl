using NCTSSoS, Test

@testset "NCTSSoS.jl" begin
    include("fastpoly_test/runtests.jl")
    include("pop.jl")
    include("sparse.jl")
    include("solver_utils.jl")
    if Sys.isapple()
        # Only test moment_problem locally due to time constraints
        include("moment_solver.jl")
    end
    include("sos_solver.jl")
    include("interface.jl")
    include("state_poly_opt.jl")
    include("Aqua.jl")
    include("Doctest.jl")
    include("ExplicitImports.jl")
end
