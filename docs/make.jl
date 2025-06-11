using Documenter
using NCTSSoS

makedocs(;
    sitename="NCTSSoS.jl",
    pages=[
        "Home" => "index.md",
        "Noncommutative Polynomial Optimization" => "ncpop.md",
        "API" => "api.md",
        "Examples" => ["Bell inequalities" => "bell.md"],
        "Performance" => ["FastPolynomials" => "performance.md"]
    ],
    modules=[NCTSSoS, NCTSSoS.FastPolynomials],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
)

deploydocs(; repo="github.com/wangjie212/NCTSSoS.git")
