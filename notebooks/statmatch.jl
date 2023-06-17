# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# +
using CSV, DataFrames

data = DataFrame(CSV.File(joinpath(@__DIR__, "..", "test", "statmatch.csv")))
# -

df = data[!, [6, 8, 9, 14, 15]]


