"""
$(SIGNATURES)

Compute a bound on the average prediction error in each base.
The bound is computed as the expected prediction error assuming that the
distribution of Z in base A (and that of Y in base B) is known, and the
prediction done with the value that maximizes the probability
"""
function bound_prediction_error(
    inst::Instance,
    norme::Int64 = 1,
    aggregate_tol::Float64 = 0.5,
)

    # Local redefinitions of parameters of  the instance
    nA = inst.nA
    nB = inst.nB
    Y = inst.Y
    Z = inst.Z

    # compute the bound in base A
    boundpredZA = 0.0
    for x = 1:length(inst.indXA)
        for mA in Y
            indwithmA = inst.indXA[x][findall(inst.Yobserv[inst.indXA[x]] .== mA)]
            nindiv = length(indwithmA)
            if nindiv == 0
                continue
            end
            for mB in Z
                indwithmB = indwithmA[findall(inst.Zobserv[indwithmA] .== mB)]
                boundpredZA +=
                    length(indwithmB) / nindiv * (1 - length(indwithmB) / nindiv) * nindiv /
                    nA
            end
        end
    end
    # @printf("Bound on average prediction error in A : %.1f %%\n", 100.*boundpredZA)

    # compute the bound in base B
    boundpredYB = 0.0
    for x = 1:length(inst.indXB)
        for mB in Z
            indwithmB = inst.indXB[x][findall(inst.Zobserv[inst.indXB[x].+nA] .== mB)] .+ nA
            nindiv = length(indwithmB)
            if nindiv == 0
                continue
            end
            for mA in Y
                indwithmA = indwithmB[findall(inst.Yobserv[indwithmB] .== mA)]
                boundpredYB +=
                    length(indwithmA) / nindiv * (1 - length(indwithmA) / nindiv) * nindiv /
                    nB
            end
        end
    end
    # @printf("Bound on average prediction error in B : %.1f %%\n", 100.*boundpredYB)

    return boundpredZA, boundpredYB

end

