"""
$(SIGNATURES)

Compute errors in the conditional distributions of a solution
"""
function compute_distrib_error!(sol::Solution, inst::Instance, empiricalZA, empiricalYB)

    nA = inst.nA
    nB = inst.nB
    Y = inst.Y
    Z = inst.Z
    nbX = length(inst.indXA)

    sol.errordistribZA = sum(
        length(inst.indXA[x][findall(inst.Yobserv[inst.indXA[x]] .== y)]) / nA *
        sum(max.(sol.estimatorZA[x, y, :] .- empiricalZA[x, y, :], 0)) for x = 1:nbX,
        y in Y
    )

    sol.errordistribYB = sum(
        length(inst.indXB[x][findall(inst.Zobserv[inst.indXB[x].+nA] .== z)]) / nB *
        sum(max.(sol.estimatorYB[x, :, z] .- empiricalYB[x, :, z], 0)) for x = 1:nbX,
        z in Z
    )

    sol.errordistribavg = (nA * sol.errordistribZA + nB * sol.errordistribYB) / (nA + nB)

end

"""
$(SIGNATURES)
"""
function compute_distrib_error_3covar(
    sol::Solution,
    inst::Instance,
    empiricalZA,
    empiricalYB,
)

    nA = inst.nA
    nB = inst.nB
    Y = inst.Y
    Z = inst.Z
    Xval = inst.Xval
    nbX = length(inst.indXA)

    Xval3 = convert(Matrix, unique(DataFrame(Xval[:, 2:3]), :auto))
    nbX3 = size(Xval3)[1]
    indX3 = Dict{Int64,Array{Int64,1}}()
    for i = 1:nbX3
        indlist = []
        for j = 1:size(inst.Xval)[1]
            if inst.Xval[j, 2:3] == Xval3[i, :]
                indlist = [indlist; j]
            end
        end
        indX3[i] = indlist
    end

    sol.errordistribZA = sum(
        max(
            sum(
                length(inst.indXA[x][findall(inst.Yobserv[inst.indXA[x]] .== y)]) / nA *
                (sol.estimatorZA[x, y, z] - empiricalZA[x, y, z]) for x in indX3[x3]
            ),
            0,
        ) for x3 = 1:nbX3, y in Y, z in Z
    )

    sol.errordistribYB = sum(
        max(
            sum(
                length(inst.indXB[x][findall(inst.Zobserv[inst.indXB[x].+nA] .== z)]) / nB * (sol.estimatorYB[x, y, z] - empiricalYB[x, y, z]) for
                x in indX3[x3]
            ),
            0,
        ) for x3 = 1:nbX3, y in Y, z in Z
    )

    sol.errordistribavg = (nA * sol.errordistribZA + nB * sol.errordistribYB) / (nA + nB)

    sol

end
