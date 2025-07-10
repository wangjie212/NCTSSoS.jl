using Aqua, NCTSSoS, Test

@testset "Aqua.jl" begin
    Aqua.test_all(NCTSSoS)
end