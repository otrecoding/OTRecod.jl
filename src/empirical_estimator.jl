"""
$(SIGNATURES)

Get an empirical estimator of the distribution of Z conditional to Y and X
on base A and reciprocally on base B obtain with a specific type of data sets

- `path`: path of the directory containing the data set
- `observed`: if nonempty, list of indices of the observed covariates; this allows to exclude some latent variables.

"""
function empirical_estimator(path, observed::Array{Int64,1} = Array{Int64,1}())

    txt_files = filter(x -> endswith(x, ".txt"), readdir(path))
    # Initialize the necessary structures by treating the first file
    # get one instance file in the direcory
    file = joinpath(path, first(txt_files))

    inst = Instance(file, 0, observed)
    nA = inst.nA
    nB = inst.nB
    Y = inst.Y
    Z = inst.Z
    indXA = inst.indXA
    indXB = inst.indXB
    nbX = length(indXA)

    # Compute the cumulative cardinality of the joint occurences of
    # (X=x,Y=y,Z=z) in the two databases
    cardA_c_mA_mB = zeros(nbX, length(Y), length(Z))
    cardB_c_mA_mB = zeros(nbX, length(Y), length(Z))
    nbIndiv = 0
    for data_file in txt_files

        # read the data file and prepare arrays
        inst = Instance(joinpath(path, data_file), 0, observed)

        # update the cumulative count
        countA, countB = empirical_distribution(inst, 0)
        cardA_c_mA_mB .+= countA
        cardB_c_mA_mB .+= countB
        nbIndiv += inst.nA
        if (nbIndiv >= 100000)
            break
        end

    end

    # Compute the empirical distribution of Z conditional to X and y in
    # base A and that of Y conditional to X and Z in base B
    # Z conditional to X and Y in base A
    empiricalZA = ones(nbX, length(Y), length(Z)) ./ length(Z)
    for x = 1:nbX, y in Y
        cardA_c_mA = sum(cardA_c_mA_mB[x, y, z] for z in Z)
        if cardA_c_mA > 0
            empiricalZA[x, y, :] .= cardA_c_mA_mB[x, y, :] ./ cardA_c_mA
        end
    end
    # Y conditional to X and Z in base B
    empiricalYB = ones(nbX, length(Y), length(Z)) ./ length(Y)
    for x = 1:nbX, z in Z
        cardB_c_mB = sum(cardB_c_mA_mB[x, y, z] for y in Y)
        if cardB_c_mB > 0
            empiricalYB[x, :, z] .= cardB_c_mA_mB[x, :, z] ./ cardB_c_mB
        end
    end

    empiricalZA, empiricalYB

end
