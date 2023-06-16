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

"""
$(SIGNATURES)

Compute the average distance between individuals of base1 with modality m1
for outcome and individuals of base2 with modality m2 for outcome

Consider only the percent_closest individuals in the computation of the
distance
"""
function avg_distance_closest(
    inst::Instance,
    base1::DataBase,
    base2::DataBase,
    outcome::DataBase,
    m1::Int,
    m2::Int,
    percent_closest::Float64,
)

    # indices of individuals of base A with given outcomes in base B (and reciprocally)
    indZinA = Dict((z, findall(inst.Zobserv[1:inst.nA] .== z)) for z in inst.Z)
    indYinB = Dict((m, findall(inst.Yobserv[inst.nA+1:end] .== y)) for y in inst.Y)

    ind1 =
        base1 == baseA ? (outcome == baseA ? inst.indY[m1] : indZinA[m1]) :
        (outcome == baseA ? indYinB[m1] : inst.indZ[m1])
    ind2 =
        base2 == baseA ? (outcome == baseA ? inst.indY[m2] : indZinA[m2]) :
        (outcome == baseA ? indYinB[m2] : inst.indZ[m2])

    # select the distance matrix depending on the base
    D =
        base1 == baseA ? (base2 == baseA ? inst.DA : inst.D) :
        (base2 == baseA ? inst.D : inst.DB)

    # swap the two sets of indices if base1=baseB and base2=baseA
    if (base1 == baseB && base2 == baseA)
        ind = ind1
        ind1 = ind2
        ind2 = ind1
    end

    # compute the average distance between the individuals in ind1 and the
    # percent_closest in ind2 and reciprocally
    avg = 0.0
    for i in ind1
        nbclose = round(Int, percent_closest * length(ind2))
        distance = sort([D[i, j] for j in ind2])
        avg += sum(distance[1:nbclose]) / nbclose / length(ind1) / 2
    end
    for j in ind2
        nbclose = round(Int, percent_closest * length(ind1))
        distance = sort([D[i, j] for i in ind1])
        avg += sum(distance[1:nbclose]) / nbclose / length(ind2) / 2
    end

    return avg
end


"""
$(SIGNATURES)

Compute prediction errors in a solution
"""
function compute_pred_error!(
    sol::Solution,
    inst::Instance,
    proba_disp::Bool = false,
    mis_disp::Bool = false,
    full_disp::Bool = false,
)

    A = 1:inst.nA
    B = 1:inst.nB
    Y = inst.Y
    Z = inst.Z
    indXA = inst.indXA
    indXB = inst.indXB
    nbX = length(indXA)

    # display the transported and real modalities
    if full_disp
        println("Modalities of base 1 individuals:")
        for i in A
            println(
                "Index: $i real value: $(inst.Zobserv[i]) transported value: $(sol.predZA[i])",
            )
        end
        # display the transported and real modalities
        println("Modalities of base 2 individuals:")
        for j in B
            println(
                "Index: $j real value: $(inst.Yobserv[inst.nA+j]) transported value: $(sol.predYB[j])",
            )
        end
    end

    # Count the number of mistakes in the transport
    #deduce the individual distributions of probability for each individual from the distributions
    probaZindivA = zeros(Float64, (inst.nA, length(Z)))
    probaYindivB = zeros(Float64, (inst.nB, length(Y)))
    for x = 1:nbX
        for i in indXA[x]
            probaZindivA[i, :] .= sol.estimatorZA[x, inst.Yobserv[i], :]
        end
        for i in indXB[x]
            probaYindivB[i, :] .= sol.estimatorYB[x, :, inst.Zobserv[i+inst.nA]]
        end
    end

    # Transport the modality that maximizes frequency
    predZA = [findmax([probaZindivA[i, z] for z in Z])[2] for i in A]
    predYB = [findmax([probaYindivB[j, y] for y in Y])[2] for j in B]

    # Base 1
    nbmisA = 0
    misA = Int64[]
    for i in A
        if predZA[i] != inst.Zobserv[i]
            nbmisA += 1
            push!(misA, i)
        end
    end

    # Base 2
    nbmisB = 0
    misB = Int64[]
    for j in B
        if predYB[j] != inst.Yobserv[inst.nA+j]
            nbmisB += 1
            push!(misB, j)
        end
    end

    if proba_disp
        if nbmisA == 0
            println("No mistake in the transport of base A")
        else
            @printf("Probability of error in base A: %.1f %%\n", 100.0 * nbmisA / inst.nA)
            if mis_disp
                println("Indices with mistakes in base A:", misA)
            end
        end

        if nbmisB == 0
            println("No mistake in the transport of base B")
        else
            @printf("Probability of error in base B: %.1f %%\n", 100.0 * nbmisB / inst.nB)
            if mis_disp
                println("Indices with mistakes in base 2:", misB)
            end
        end
    end

    sol.errorpredZA = nbmisA / inst.nA
    sol.errorpredYB = nbmisB / inst.nB
    sol.errorpredavg =
        (inst.nA * sol.errorpredZA + inst.nB * sol.errorpredYB) / (inst.nA + inst.nB)

end

"""
$(SIGNATURES)

Compute errors in the conditional distributions of a solution
"""
function compute_distrib_error(inst::Instance, sol::Solution, empiricalZA, empiricalYB)

    nA = inst.nA
    nB = inst.nB
    Y = copy(inst.Y)
    Z = copy(inst.Z)
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

    sol

end

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

"""
$(SIGNATURES)

Compute the cost between pairs of outcomes as the average distance between
covariations of individuals with these outcomes, but considering only the
percent closest neighbors
"""
function average_distance_to_closest(inst::Instance, percent_closest::Float64)

    # Redefine A and B for the model
    A = 1:inst.nA
    B = 1:inst.nB
    Y = inst.Y
    Z = inst.Z
    indY = inst.indY
    indZ = inst.indZ

    # Compute average distances as described in the above
    Davg = zeros(Float64, (length(Y), length(Z)))
    DindivA = zeros(Float64, (inst.nA, length(Z)))
    DindivB = zeros(Float64, (inst.nB, length(Y)))

    for y in Y, i in indY[y], z in Z

        nbclose = max(round(Int, percent_closest * length(indZ[z])), 1)
        distance = sort([inst.D[i, j] for j in indZ[z]])
        DindivA[i, z] = sum(distance[1:nbclose]) / nbclose
        Davg[y, z] += sum(distance[1:nbclose]) / nbclose / length(indY[y]) / 2.0

    end

    for z in Z, j in indZ[z], y in Y

        nbclose = max(round(Int, percent_closest * length(indY[y])), 1)
        distance = sort([inst.D[i, j] for i in indY[y]])
        DindivB[j, y] = sum(distance[1:nbclose]) / nbclose
        Davg[y, z] += sum(distance[1:nbclose]) / nbclose / length(indZ[z]) / 2.0

    end

    Davg, DindivA, DindivB

end

"""
$(SIGNATURES)

Compute prediction errors in a solution
"""
function compute_pred_error(
    inst::Instance,
    sol::Solution,
    proba_disp::Bool = false,
    mis_disp::Bool = false,
    full_disp::Bool = false,
)

    A = 1:inst.nA
    B = 1:inst.nB
    Y = inst.Y
    Z = inst.Z
    indXA = inst.indXA
    indXB = inst.indXB
    nbX = length(indXA)

    # display the transported and real modalities
    if full_disp
        println("Modalities of base 1 individuals:")
        for i in A
            println(
                "Index: $i real value: $(inst.Zobserv[i]) transported value: $(sol.predZA[i])",
            )
        end
        # display the transported and real modalities
        println("Modalities of base 2 individuals:")
        for j in B
            println(
                "Index: $j real value: $(inst.Yobserv[inst.nA+j]) transported value: $(sol.predYB[j])",
            )
        end
    end

    # Count the number of mistakes in the transport
    #deduce the individual distributions of probability for each individual from the distributions
    probaZindivA = Array{Float64,2}(undef, inst.nA, length(Z))
    probaYindivB = Array{Float64,2}(undef, inst.nB, length(Y))
    for x = 1:nbX
        for i in indXA[x]
            probaZindivA[i, :] = sol.estimatorZA[x, inst.Yobserv[i], :]
        end
        for i in indXB[x]
            probaYindivB[i, :] = sol.estimatorYB[x, :, inst.Zobserv[i+inst.nA]]
        end
    end

    # Transport the modality that maximizes frequency
    predZA = [findmax([probaZindivA[i, z] for z in Z])[2] for i in A]
    predYB = [findmax([probaYindivB[j, y] for y in Y])[2] for j in B]

    # Base 1
    nbmisA = 0
    misA = []
    for i in A
        if predZA[i] != inst.Zobserv[i]
            nbmisA += 1
            push!(misA, i)
        end
    end

    # Base 2
    nbmisB = 0
    misB = []
    for j in B
        if predYB[j] != inst.Yobserv[inst.nA+j]
            nbmisB += 1
            push!(misB, j)
        end
    end

    if proba_disp
        if nbmisA == 0
            println("No mistake in the transport of base A")
        else
            @printf("Probability of error in base A: %.1f %%\n", 100.0 * nbmisA / inst.nA)
            if mis_disp
                println("Indices with mistakes in base A:", misA)
            end
        end

        if nbmisB == 0
            println("No mistake in the transport of base B")
        else
            @printf("Probability of error in base B: %.1f %%\n", 100.0 * nbmisB / inst.nB)
            if mis_disp
                println("Indices with mistakes in base 2:", misB)
            end
        end
    end

    sol.errorpredZA, sol.errorpredYB = nbmisA / inst.nA, nbmisB / inst.nB
    sol.errorpredavg =
        (inst.nA * sol.errorpredZA + inst.nB * sol.errorpredYB) / (inst.nA + inst.nB)

    sol

end
