using Documenter
using NCTSSoS
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:authoryear)

makedocs(;
    sitename="NCTSSoS.jl",
    pages=[
        "Home" => "index.md",
        "Quick Start" => "quick_start.md",
        "Backgrounds" => "backgrounds.md",
        "Manual" => ["Polynomials"=> "manual/polynomials.md","Sparsities" => "manual/sparsities.md", "SDP Relaxation"=> "manual/sdp_relaxation.md"],
        "Examples" => ["Bell inequalities" => "examples/bell.md", "Certifying Ground State" => "examples/cert_ground_state.md",
            "Noncommutative Polynomial Optimization" => "examples/ncpop.md",
            ],
        "Optimizers" => "optimizers.md",
        "References" => "reference.md",
        "APIs" => ["User interface" => "apis/interface.md", "Polynomials"=> "apis/polynomials.md", "Sparsities"=> "apis/sparsities.md", "SDP Relaxation"=> "apis/relaxations.md"]],
        plugins=[bib],
        modules=[NCTSSoS, NCTSSoS.FastPolynomials],
        format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true", assets=String["assets/citations.css"]),
)

deploydocs(; repo="github.com/wangjie212/NCTSSoS.git", devbranch="main")
