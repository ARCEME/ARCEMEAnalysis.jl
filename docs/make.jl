using ARCEMEAnalysis
using Documenter

DocMeta.setdocmeta!(ARCEMEAnalysis, :DocTestSetup, :(using ARCEMEAnalysis); recursive=true)

makedocs(;
    modules=[ARCEMEAnalysis],
    authors="Melanie Weynants <mweynants@bgc-jena.mpg.de>, Fabian Gans <fgans@bgc-jena.mpg.de>",
    sitename="ARCEMEAnalysis.jl",
    format=Documenter.HTML(;
        canonical="https://ARCEME.github.io/ARCEMEAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ARCEME/ARCEMEAnalysis.jl",
    devbranch="main",
)
