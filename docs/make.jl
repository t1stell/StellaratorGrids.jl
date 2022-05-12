using StellaratorGrids
using Documenter

DocMeta.setdocmeta!(StellaratorGrids, :DocTestSetup, :(using StellaratorGrids); recursive=true)

makedocs(;
    modules=[StellaratorGrids],
    authors="Benjamin Faber <bfaber@wisc.edu> and contributors",
    repo="https://gitlab.com/bfaber/StellaratorGrids.jl/blob/{commit}{path}#{line}",
    sitename="StellaratorGrids.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bfaber.gitlab.io/StellaratorGrids.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
