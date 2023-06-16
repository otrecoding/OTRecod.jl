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

