using TreeContractor
using Documenter

DocMeta.setdocmeta!(TreeContractor, :DocTestSetup, :(using TreeContractor); recursive=true)

makedocs(;
    modules=[TreeContractor],
    authors="nzy1997",
    sitename="TreeContractor.jl",
    format=Documenter.HTML(;
        canonical="https://nzy1997.github.io/TreeContractor.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nzy1997/TreeContractor.jl",
    devbranch="main",
)
