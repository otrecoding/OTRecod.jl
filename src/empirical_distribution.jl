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

