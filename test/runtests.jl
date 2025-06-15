using NCTSSoS, Test

@testset "NCTSSoS.jl" begin
    include("pop.jl")
    include("sparse.jl")
    # TODO:
    include("solver_utils.jl")
    include("moment_solver.jl")
    include("sos_solver.jl")
    include("interface.jl")
    include("state_moment_solver.jl")
end
