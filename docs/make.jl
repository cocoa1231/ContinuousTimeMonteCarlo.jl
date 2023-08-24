using ContinuousTimeMonteCarlo
using Documenter

DocMeta.setdocmeta!(ContinuousTimeMonteCarlo, :DocTestSetup, :(using ContinuousTimeMonteCarlo); recursive=true)

makedocs(;
    modules=[ContinuousTimeMonteCarlo],
    authors="Cocoa K.",
    repo="https://github.com/cocoa1231/ContinuousTimeMonteCarlo.jl/blob/{commit}{path}#{line}",
    sitename="ContinuousTimeMonteCarlo.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cocoa1231.github.io/ContinuousTimeMonteCarlo.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cocoa1231/ContinuousTimeMonteCarlo.jl",
    devbranch="master",
)
