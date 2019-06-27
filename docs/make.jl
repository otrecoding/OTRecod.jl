using Documenter
using OTRecod

makedocs(
    sitename = "OTRecod",
    format = Documenter.HTML(),
    modules = [OTRecod]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
