using Robin2011_RepPackage
using Documenter

DocMeta.setdocmeta!(Robin2011_RepPackage, :DocTestSetup, :(using Robin2011_RepPackage); recursive=true)

makedocs(;
    modules=[Robin2011_RepPackage],
    authors="Bo Jacobs Strom <bojs.13@gmail.com> and contributors",
    sitename="Robin2011_RepPackage.jl",
    format=Documenter.HTML(;
        canonical="https://bo-js.github.io/Robin2011_RepPackage.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/bo-js/Robin2011_RepPackage.jl",
    devbranch="main",
)
