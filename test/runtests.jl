using NCTSSoS, Test

@testset "FastPolynomials.jl" begin
    include("fastpoly_test/runtests.jl")
end

@testset "pop.jl" begin
    include("pop.jl")
end

@testset "sparse.jl" begin
    include("sparse.jl")
end

@testset "solver_utils.jl" begin
    include("solver_utils.jl")
end

@testset "moment_solver.jl" begin
    if Sys.isapple()
        # Only test moment_problem locally due to time constraints
        include("moment_solver.jl")
    end
end

@testset "sos_solver.jl" begin
    include("sos_solver.jl")
end

@testset "interface.jl" begin
    include("interface.jl")
end

@testset "state_poly_opt.jl" begin
    include("state_poly_opt.jl")
end

@testset "trace_poly_opt.jl" begin
    include("trace_poly_opt.jl")
end

@testset "Aqua.jl" begin
    include("Aqua.jl")
end

@testset "Doctest.jl" begin
    include("Doctest.jl")
end

@testset "ExplicitImports.jl" begin
    include("ExplicitImports.jl")
end