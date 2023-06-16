@enum DataBase baseA baseB

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

"""
$(SIGNATURES)

Return the empirical cardinality of the joint occurrences of (C=x,Y=mA,Z=mB)
in both bases
"""
function empirical_distribution(
    inst::Instance,
    norme::Int64 = 0,
    aggregate_tol::Float64 = 0.5,
)

    # Local redefinitions of parameters of  the instance
    Y = inst.Y
    Z = inst.Z

    # aggregate the individuals per covariate
    nbX = length(inst.indXA)

    # count the cardinality of occurrence of each triplet (x,mA,mB) in both bases
    cardA_c_mA_mB = zeros(nbX, length(Y), length(Z))
    for x = 1:nbX
        for i in inst.indXA[x]
            cardA_c_mA_mB[x, inst.Yobserv[i], inst.Zobserv[i]] += 1
        end
    end
    cardB_c_mA_mB = zeros(nbX, length(Y), length(Z))
    for x = 1:nbX
        for i in inst.indXB[x]
            cardB_c_mA_mB[x, inst.Yobserv[i+inst.nA], inst.Zobserv[i+inst.nA]] += 1
        end
    end

    return cardA_c_mA_mB, cardB_c_mA_mB
end


"""
$(SIGNATURES)

Display information about the distance between the modalities
"""
function disp_inst_info(inst::Instance)

    # local definitions
    nA = inst.nA
    nB = inst.nB
    A = 1:inst.nA
    B = 1:inst.nB
    Y = inst.Y
    Z = inst.Z
    indY = inst.indY
    indZ = inst.indZ


    println("\n#################################################################")
    println("INFORMATION ABOUT THE INSTANCE")
    println("#################################################################\n")

    # return indicators about the original density of the modalities
    print("Average distance between objects of base 1: ")
    @printf("%.2f\n", 10 * sum([DA[i, j] for i in A, j in A]) / (nA^2))
    print("Average distance between objects of base 2: ")
    @printf("%.2f\n", 10 * sum([DB[i, j] for i in B, j in B]) / (nB^2))
    print("Crossed average distance between objects of base 1 and 2: ")
    @printf("%.2f\n", 10 * sum([inst.D[i, j] for i in A, j in B]) / (nA * nB))

    # restrict the average distance to the 10% closest individuals
    println("\nAverage distance between objects per modality")
    println("Modalities of Y, individuals of A:")
    percent_closest = 0.1
    for y1 in Y
        for y2 in Y
            if y1 > y2
                continue
            end
            avg = avg_distance_closest(inst, baseA, baseA, baseA, y1, y2, 1.0)
            @printf("\tModalities %d and %d : %.2f\n", y1, y2, 10 * avg)
            avg = avg_distance_closest(inst, baseA, baseA, baseA, y1, y2, percent_closest)
            @printf(
                "\t\trestricted to the %.1f %% closest: %.2f\n",
                100 * percent_closest,
                10 * avg
            )
        end
    end
    println("\nModalities of Z, individuals of B:")
    for z1 in Z
        for z2 in Z
            if z1 > z2
                continue
            end
            avg = avg_distance_closest(inst, baseB, baseB, baseB, z1, z2, 1.0)
            @printf("\tModalities %d and %d : %.2f\n", z1, z2, 10 * avg)
            avg = avg_distance_closest(inst, baseB, baseB, baseB, z1, z2, percent_closest)
            @printf(
                "\t\trestricted to the %.1f %% closest: %.2f\n",
                100 * percent_closest,
                10 * avg
            )
        end
    end
    println("\nModalities of Y, crossed bases:")
    for y1 in Y
        for y2 in Y
            avg = avg_distance_closest(inst, baseA, baseB, baseA, y1, y2, 1.0)
            @printf("\tModalities %d and %d : %.2f\n", y1, y2, 10 * avg)
            avg = avg_distance_closest(inst, baseA, baseB, baseA, y1, y2, percent_closest)
            @printf(
                "\t\trestricted to the %.1f %% closest: %.2f\n",
                100 * percent_closest,
                10 * avg
            )
        end
    end
    println("\nModalities of Z, crossed bases:")
    for z1 in Z
        for z2 in Z
            avg = avg_distance_closest(inst, baseA, baseB, baseB, z1, z2, 1.0)
            @printf("\tModalities %d and %d : %.2f\n", z1, z2, 10 * avg)
            avg = avg_distance_closest(inst, baseA, baseB, baseB, z1, z2, percent_closest)
            @printf(
                "\t\trestricted to the %.1f%% closest: %.2f\n",
                100.0 * percent_closest,
                10 * avg
            )
        end
    end
end

