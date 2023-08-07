push!(LOAD_PATH, "../src/")
using AngularCorrelation

using Documenter
makedocs(
    sitename="AngularCorrelation.jl",
    modules=[AngularCorrelation],
    pages=[
        "AngularCorrelation.jl" => "index.md"
    ])
deploydocs(;
    repo="https://github.com/op3/AngularCorrelation.jl"
)
