"""
$(SIGNATURES)
"""
function aggregate_per_covar_mixed(
    inst::Instance,
    norme::Int64 = 1,
    aggregate_tol::Float64 = 0.5,
)

    A = 1:inst.nA
    B = 1:inst.nB

    # initialization of the structures
    nbX = 0
    indXA = Dict{Int64,Array{Int64}}()
    indXB = Dict{Int64,Array{Int64}}()
    notaggA = [i for i in A]
    notaggB = [i for i in B]

    # aggregate until every individual in base A is aggregated
    while !isempty(notaggA)
        nbX += 1
        ind = notaggA[1]
        isinset = inst.DA[ind, notaggA] .< aggregate_tol
        indXA[nbX] = notaggA[isinset]
        deleteat!(notaggA, isinset)
        isinset = inst.D[ind, notaggB] .< aggregate_tol
        indXB[nbX] = notaggB[isinset]
        deleteat!(notaggB, isinset)
    end

    # complete the aggregation with the individuals of base B that are not aggregated yet
    while !isempty(notaggB)
        nbX += 1
        ind = notaggB[1]
        isinset = inst.DB[ind, notaggB] .< aggregate_tol
        indXB[nbX] = notaggB[isinset]
        indXA[nbX] = []
        deleteat!(notaggB, isinset)
    end

    return indXA, indXB
end

