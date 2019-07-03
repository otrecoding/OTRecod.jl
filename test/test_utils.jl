using CSV, Distances

@enum IndivFromGroup sequential optimal

@info "... compute the empirical distributions of outcomes"

path = joinpath("data","tab.csv")

df = CSV.read(path)

maxrelax = 0.0
percent_closest = 0.2

for method in [group, joint]

    for (norme,distance) in enumerate([Cityblock(), Euclidean(), Hamming()])

        inst  = Instance(df, distance)

        nA    = inst.nA
        nB    = inst.nB
        Y     = copy(inst.Y)
        Z     = copy(inst.Z)
        indXA = copy(inst.indXA)
        indXB = copy(inst.indXB)
        nbX   = length(indXA)

        @debug """ Compute the cumulative cardinality of 
                   the joint occurences of 
                   (X=x,Y=y,Z=z) in the two databases """

        cardA_c_mA_mB = zeros(nbX,length(Y),length(Z))
        cardB_c_mA_mB = zeros(nbX,length(Y),length(Z))
        nbIndiv = 0

        # Compute the empirical distribution of Z conditional 
        # to X and y in base A and that of Y conditional to X 
        # and Z in base B Z conditional to X and Y in base A

        empiricalZA = 1/length(Z)*ones(nbX,length(Y),length(Z))
        for x in 1:nbX, y in Y
            cardA_c_mA = sum(cardA_c_mA_mB[x,y,z] for z in Z)
            if cardA_c_mA > 0
                empiricalZA[x,y,:] = cardA_c_mA_mB[x,y,:]/cardA_c_mA
            end
        end
        # Y conditional to X and Z in base B
        empiricalYB = 1/length(Y)*ones(nbX,length(Y),length(Z))
        for x in 1:nbX, z in Z
            cardB_c_mB = sum(cardB_c_mA_mB[x,y,z] for y in Y)
            if cardB_c_mB > 0
                empiricalYB[x,:,z] = cardB_c_mA_mB[x,:,z]/cardB_c_mB
            end
        end

        # compute a bound on the average prediction error in each base
        # println("... compute bounds on the average prediction errors\n")
        # errorboundA, errorboundB = compute_average_error_bound(path, norme)

        nA    = inst.nA 
        nB    = inst.nB 
        Y     = copy(inst.Y)
        Z     = copy(inst.Z)
        indXA = copy(inst.indXA)
        indXB = copy(inst.indXB)
        nbX   = length(indXA)

        if method == group
           indiv_method = maxrelax > 0.0 ? optimal : sequential
           sol = OT_group(inst,percent_closest,maxrelax,norme,indiv_method)
           lambda_reg = 0.0;
        elseif method == joint
           sol = OT_joint(inst, maxrelax, lambda_reg, percent_closest)
        end

        sol = OTRecod.compute_pred_error(inst, sol, false);
        sol = OTRecod.compute_distrib_error(inst, sol, empiricalZA, 
                                            empiricalYB)

        @printf("%-12s , %-6s , %-5d , %-5.2f , %-5.2f , %-6.3f , %-6.3f , %-8.3f , %-7.3f , %-7.3f , %-9.3f , %-6.2f\n", data_file, method, norme, maxrelax, lambda_reg, sol.errorpredZA, sol.errorpredYB, sol.errorpredavg, sol.errordistribZA, sol.errordistribYB, sol.errordistribavg, sol.tsolve)

    end

end
