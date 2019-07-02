using Documenter
using OTRecod

makedocs(
    sitename = "OTRecod",
    authors  = "Jeremy Omer",
    format   = Documenter.HTML(),
    modules  = [OTRecod],
    pages    = [ "Documentation" => "index.md",
                 "Group"         => "OT_group.md",
                 "Joint"         => "OT_joint.md",
                 "Logging"       => "PrintLog.md",
                 "Plotting"      => "plot_functions.md",
                 "Utilities"     => "utils.md"]
)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/otrecoding/OTRecod.jl.git",
 )
