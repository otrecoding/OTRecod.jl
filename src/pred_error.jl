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
