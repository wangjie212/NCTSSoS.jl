using Test, Documenter, NCTSSoS

@testset "DocTest" begin
    DocMeta.setdocmeta!(NCTSSoS, :DocTestSetup, :(using NCTSSoS, Clarabel, COSMO); recursive=true)
    doctest(NCTSSoS; manual=false, fix=false)
end