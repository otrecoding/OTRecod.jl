using Documenter
using OTRecoding

makedocs(
    sitename = "OTRecoding",
    format = Documenter.HTML(),
    modules = [OTRecoding]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
