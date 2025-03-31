using NCTSSoS, Test

@testset "NCTSSoS.jl" begin
    include("pop.jl")
    include("sparse.jl")
    include("solver_utils.jl")
    include("moment_solver.jl")
    include("sos_solver.jl")
    include("interface.jl")
end
