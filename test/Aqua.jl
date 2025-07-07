using Aqua, NCTSSoS, Test

@testset "Aqua.jl" begin
    Aqua.test_all(NCTSSoS;
        stale_deps=(ignore = [:JET, :DispatchDoctor],)
    )
end