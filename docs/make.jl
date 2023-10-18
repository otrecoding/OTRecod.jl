push!(LOAD_PATH, "../src/")

using Documenter
using OTRecod

makedocs(
    sitename = "OTRecod",
    authors = "Jeremy Omer",
    format = Documenter.HTML(),
    modules = [OTRecod],
    pages = [
        "Documentation" => "index.md",
        "Group" => "ot_group.md",
        "Joint" => "ot_joint.md",
        "Logging" => "PrintLog.md",
        "Utilities" => "utils.md",
    ],
)

deploydocs(
    repo = "github.com/otrecoding/OTRecod.jl.git",
)
