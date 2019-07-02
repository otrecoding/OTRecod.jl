using CSV

@info "... compute the empirical distributions of outcomes"

path = joinpath("data","tab.csv")

df = DataFrame(CSV.File("test/data/tab.csv"))

norme = 1
#=

files = readdir(path);
dir_name = split(path,"/")[end];

## Initialize the necessary structures by treating the first file
file = "file.txt";
# get one instance file in the direcory
for data_file in files
    # continue if not a data file
    if (data_file[end-3:end]==".txt")
        file = string(path,"/",data_file);
        break;
    end
end
inst  = Instance(file, norme);
nA    = inst.nA
nB    = inst.nB
Y     = copy(inst.Y)
Z     = copy(inst.Z)
indXA = copy(inst.indXA)
indXB = copy(inst.indXB)
nbX   = length(indXA)

## Compute the cumulative cardinality of the joint occurences of (X=x,Y=y,Z=z) in the two databases
cardA_c_mA_mB = zeros(nbX,length(Y),length(Z));
cardB_c_mA_mB = zeros(nbX,length(Y),length(Z));
nbIndiv = 0;
for data_file in files
    # continue if not a data file
    if !(data_file[end-3:end]==".txt") continue end

    # read the data file and prepare arrays
    inst = Instance(string(path,"/",data_file), norme);

    # update the cumulative count
    countA, countB = empirical_distribution(inst, norme);
    cardA_c_mA_mB = cardA_c_mA_mB .+ countA;
    cardB_c_mA_mB = cardB_c_mA_mB .+ countB;
    nbIndiv += inst.nA;
    if (nbIndiv >= 100000) break end
end

## Compute the empirical distribution of Z conditional to X and y in base A and that of Y conditional to X and Z in base B
# Z conditional to X and Y in base A
empiricalZA = 1/length(Z)*ones(nbX,length(Y),length(Z));
for x in 1:nbX
    for y in Y
        cardA_c_mA = sum(cardA_c_mA_mB[x,y,z] for z in Z);
        if cardA_c_mA > 0
            empiricalZA[x,y,:] = cardA_c_mA_mB[x,y,:]/cardA_c_mA;
        end
    end
end
# Y conditional to X and Z in base B
empiricalYB = 1/length(Y)*ones(nbX,length(Y),length(Z));
for x in 1:nbX
    for z in Z
        cardB_c_mB = sum(cardB_c_mA_mB[x,y,z] for y in Y);
        if cardB_c_mB > 0
            empiricalYB[x,:,z] = cardB_c_mA_mB[x,:,z]/cardB_c_mB;
        end
    end
end


    # compute a bound on the average prediction error in each base
    # println("... compute bounds on the average prediction errors\n");
    # errorboundA, errorboundB = compute_average_error_bound(path, norme);

    # solve the instances corresponding to each file
    files = readdir(path)
    nbruns = 0;
    println("Maxrelax= ", maxrelax);
    for data_file in files
        # stop if the requested number of runs has been performed
        if ((nbfiles > 0) & (nbruns >= nbfiles)) break end
        # continue if not a data file
        if !(data_file[end-3:end]==".txt") continue end

        # Reading the data file and preparing arrays
        inst = Instance(string(path,"/",data_file), norme);
        nA = inst.nA ;  nB = inst.nB ; Y = copy(inst.Y); Z = copy(inst.Z);
        indXA = copy(inst.indXA); indXB = copy(inst.indXB);
        nbX = length(indXA);

        println("\n########## File : ", string(path,"/",data_file), " ##########");
        if method == group
            indiv_method = maxrelax > 0.0 ? optimal : sequential;
            sol = OT_group(inst,percent_closest,maxrelax,norme,indiv_method);
            lambda_reg = 0.0;
        elseif method == joint
            sol = OT_joint(inst, maxrelax, lambda_reg, percent_closest);
        end

        sol = compute_pred_error(inst, sol, false);
        sol = compute_distrib_error(inst, sol, empiricalZA, empiricalYB);

        outfile = open(outname, "a");
        @printf(outfile, "%-12s , %-6s , %-5d , %-5.2f , %-5.2f , %-6.3f , %-6.3f , %-8.3f , %-7.3f , %-7.3f , %-9.3f , %-6.2f\n", data_file, method, norme, maxrelax, lambda_reg, sol.errorpredZA, sol.errorpredYB, sol.errorpredavg, sol.errordistribZA, sol.errordistribYB, sol.errordistribavg, sol.tsolve);
        close(outfile);

        nbruns += 1
    end
end

=#
