using CSV, Distances, MLDataUtils
using CategoricalArrays

@testset "StatMatch" begin

    df = CSV.read("statmatch.csv", copycols = true)
    df = df[[6, 8, 9, 14, 15]]

    covariables = [:age, :hsize, :sex]
    outcomes = [:Y, :Z]

    # rename columns
    names!(df, vcat(covariables, outcomes))

    categorical!(df, compress = true)
    # clean data

    maxrelax = 0.0
    percent_closest = 0.2
    lambda_reg = 0.0

    for method in [:group, :joint]

        for distance in [Cityblock(), Euclidean(), Hamming()]

            instance = Instance(df, covariables, outcomes, distance)

            #=
                        nA    = instance.nA
                        nB    = instance.nB
                        Y     = instance.Y
                        Z     = instance.Z
                        indXA = instance.indXA
                        indXB = instance.indXB
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
                                empiricalZA[x,y,:] .= cardA_c_mA_mB[x,y,:] ./ cardA_c_mA
                            end
                        end
                        # Y conditional to X and Z in base B
                        empiricalYB = 1/length(Y)*ones(nbX,length(Y),length(Z))
                        for x in 1:nbX, z in Z
                            cardB_c_mB = sum(cardB_c_mA_mB[x,y,z] for y in Y)
                            if cardB_c_mB > 0
                                empiricalYB[x,:,z] .= cardB_c_mA_mB[x,:,z] ./ cardB_c_mB
                            end
                        end

                        # compute a bound on the average prediction error in each base
                        # println("... compute bounds on the average prediction errors\n")
                        # errorboundA, errorboundB = compute_average_error_bound(path, norme)

                        if method == :group

                           indiv_method = maxrelax > 0.0 ? :optimal : :sequential

                           sol = ot_group(instance, percent_closest, 
                                          maxrelax, indiv_method)

                        elseif method == :joint

                           sol = ot_joint(instance, maxrelax, 
                                          lambda_reg, percent_closest)

                        end

                        OTRecod.compute_pred_error!( sol, instance, false)
                        OTRecod.compute_distrib_error!(sol, instance, 
                                                       empiricalZA, 
                                                       empiricalYB)

            =#

            @show method, norme
            @show sol
            @test true

        end

    end


end
