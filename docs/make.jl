using Documenter
using NCTSSoS
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:authoryear)

makedocs(;
    sitename="NCTSSoS.jl",
    pages=[
        "Home" => "index.md",
        "Noncommutative Polynomial Optimization" => "ncpop.md",
        "API" => "api.md",
        "Examples" => ["Bell inequalities" => "bell.md", "Broyden Banded Function" => "broyden.md", "Certifying Ground State" => "cert_ground_state.md"],
        "References" => "reference.md"
    ],
    plugins=[bib],
    modules=[NCTSSoS, NCTSSoS.FastPolynomials],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true", assets=String["assets/citations.css"]),
)

deploydocs(; repo="github.com/wangjie212/NCTSSoS.git", devbranch="main")
