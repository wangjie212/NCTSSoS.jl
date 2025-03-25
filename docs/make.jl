using Documenter
using NCTSSoS

makedocs(;
    sitename="NCTSSoS.jl",
    pages=[
        "Home" => "index.md",
        "API" => "api.md"
    ],
    modules=[NCTSSOS],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
)

deploydocs(; repo="github.com/wangjie212/NCTSSoS.git")
