# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.1
#   kernelspec:
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# +
using Parameters
using Distributions

digitize(x, bins) = searchsortedlast.(Ref(bins), x)

function to_categorical(x)
    res = sort(unique(x)) .== reshape(x, (1, size(x)...))
    return res[2:end,:]
end
# -

@with_kw struct Data
    nA :: Int = 1000
    nB :: Int = 1000
    mA :: Vector{Float64} = [0, 0, 0]
    mB :: Vector{Float64} = [0, 0, 0]
    covA :: Matrix{Float64} = [1.0 0.2 0.2; 0.2 1.0 0.2; 0.2 0.2 1.0]
    covB :: Matrix{Float64} = [1.0 0.2 0.2; 0.2 1.0 0.2; 0.2 0.2 1.0]
    px1c :: Vector{Float64} = [0.5, 0.5]
    px2c :: Vector{Float64} = [0.333, 0.334, 0.333]
    px3c :: Vector{Float64} = [0.25, 0.25, 0.25, 0.25]
    p :: Float64 =0.6
    aA :: Vector{Float64} = [1, 1, 1, 1, 1, 1]
    aB :: Vector{Float64} = [1, 1, 1, 1, 1, 1]
    eps :: Float64 = 0
end

data_parameters = Data()

# +
d = MvNormal(data_parameters.mA, data_parameters.covA)    
X = rand(d, data_parameters.nA)

px1cc = cumsum(data_parameters.px1c)[1:end-1]
px2cc = cumsum(data_parameters.px2c)[1:end-1]
px3cc = cumsum(data_parameters.px3c)[1:end-1]

qx1c = quantile(Normal(0.0, 1.0), px1cc)
qx2c = quantile(Normal(0.0, 1.0), px2cc)
qx3c = quantile(Normal(0.0, 1.0), px3cc)

bins11 = vcat(minimum(Xglob1[1,:]) - 100, qx1c, maximum(X[1,:]) + 100)
bins12 = vcat(minimum(Xglob1[2,:]) - 100, qx2c, maximum(X[2,:]) + 100)
bins13 = vcat(minimum(Xglob1[3,:]) - 100, qx3c, maximum(X[3,:]) + 100)

X = rand(d, data_parameters.nA)

X1 = digitize(X[1,:], bins11)
X2 = digitize(X[2,:], bins12)
X3 = digitize(X[3,:], bins13)

X1c = to_categorical(X1)
X2c = to_categorical(X2)
X3c = to_categorical(X3)
# -

Y1 = vcat(X1c, X2c, X3c)' * data_parameters.aA

Y1 = X_dumm1




