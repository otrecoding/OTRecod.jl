using Documenter
using OTRecod

makedocs(
    sitename = "OTRecod",
    authors  = "Jeremy Omer",
    format   = Documenter.HTML(),
    modules  = [OTRecod]
)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/otrecoding/OTRecod.jl.git",
 )
