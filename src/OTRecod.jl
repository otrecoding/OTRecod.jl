module OTRecod

using LaTeXStrings
using Statistics
using JuMP, Cbc, Clp
using Printf
using StatsPlots

export run_directory

include("utils.jl")
include("OT_group.jl")
include("OT_joint.jl")
include("plot_functions.jl")


"""
    compute_average_error_bound(path; norme=1)

Compute a lower bound on the best average prediction error that one can
obtain with a specific type of data sets
path: path of the directory containing the data set
"""
function compute_average_error_bound(path, norme::Int64=1)

    files = readdir(path);
    dir_name = split(path,"/")[end];
    errorboundA = 0.0;
    errorboundB = 0.0;
    nbfiles = 0;

    for data_file in files
        # continue if not a data file
        if !(data_file[end-3:end]==".txt") continue end

        # Reading the data file and preparing arrays
        inst = Instance(string(path,"/",data_file), norme);

        # compute the bound and update the cumulative value
        boundA,boundB = bound_prediction_error(inst,norme);
        errorboundA += boundA;
        errorboundB += boundB;
        nbfiles += 1;
    end
    errorboundA = errorboundA/nbfiles;
    errorboundB = errorboundB/nbfiles;

    @printf("Bound on average prediction error in A : %.1f %%\n", 100.0 * errorboundA);
    @printf("Bound on average prediction error in B : %.1f %%\n", 100.0 * errorboundB);

    return errorboundA, errorboundB;
end

"""
    empirical_estimator(path; norme=1)

Get an empirical estimator of the distribution of Z conditional to Y and X
on base A and reciprocally on base B
obtain with a specific type of data sets
path: path of the directory containing the data set
"""
function empirical_estimator(path, norme::Int64=1)

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
    inst = Instance(file, norme);
    nA = inst.nA; nB = inst.nB; Y = copy(inst.Y); Z = copy(inst.Z);
    indXA = copy(inst.indXA); indXB = copy(inst.indXB);
    nbX = length(indXA);

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

    return empiricalZA, empiricalYB;

end

"""
    run_directory(path, method; outname="result.out", 
                                maxrelax=0.0, 
                                lambda_reg=0.0, 
                                nbfiles=0, 
                                norme=0, 
                                percent_closest=0.2)

 Run one given method on a given number of data files of a given directory
 The data files must be the only files with extension ".txt" in the directory

 - `path`   : name of the directory
 - `method` : `:group` or `:joint`
 - `nbfiles`: number of files considered, 0 if all the data files are tested
 - `norme`  : 1 or 2, norm used for distances in the space of covariates

 (see run_all_methods for the description of other parameters)
"""
function run_directory(path, method, outname::String="result.out", 
                       maxrelax::Float64=0.0, lambda_reg::Float64=0.0, 
                       nbfiles::Int64=0, norme::Int64=0, 
                       percent_closest::Float64=0.2)

    println("\n#################################################################")
    println("RUN ONE METHOD ON THE FILES OF ONE DIRECTORY ")
    println("\tMethod: ", method);
    println("\tDirectory: ", path);
    println("\tOutput file: ", outname);
    if (nbfiles > 0) println("\tTest only ", nbfiles, " files");
    else println("\tTest all the files of the directory");end
    if (method == :joint) println("\tRegularization parameter: ", lambda_reg);end
    println("\n#################################################################\n")

    # initialize the output file
    outfile = open(outname, "w");
    @printf(outfile, "%-12s , %-6s , %-5s , %-5s , %-5s , %-6s , %-6s , %-8s , %-7s , %-7s , %-9s , %-6s\n","filename", "method", "norme", "relax", "regul", "ePredA", "ePredB", "ePredavg", "eDistrA", "eDistrB","eDistravg","cpu");
    close(outfile);

    # compute the true empirical conditional distributions for comparison with results
    println("... compute the empirical distributions of outcomes\n");
    empiricalZA,empiricalYB = empirical_estimator(path, norme);

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
        if method == :group
            indiv_method = maxrelax > 0.0 ? optimal : sequential;
            sol = OT_group(inst,percent_closest,maxrelax,norme,indiv_method);
            #PN lambda_reg = 0.0;
        elseif method == :joint
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



"""
    run_benchmark(path, method, maxrelax=0.0, lambda_reg=0.0, norme=0, 
                  percent_closest=0.2)

 Run one method on the complete benchmark
 path: path of the directory including the benchmark
"""
function run_benchmark(path, method, maxrelax::Float64=0.0, 
                       lambda_reg::Float64=0.0, norme::Int64=0,
                       percent_closest::Float64=0.2)

    println("\n#################################################################")
    println("RUN ONE METHOD ON THE COMPLETE BENCHMARK ")
    println("\tMethod: ", method);
    if (method == :joint) println("\tRegularization parameter: ", lambda_reg);end
    println("\n#################################################################\n")

    dirlist = readdir(path);
    restart = true
    for dir in dirlist
        datasetpath = string(path,"/",dir);
        println(datasetpath);
        if !isdir(datasetpath) continue end
        if (dir != "Sn-250") && (dir != "Sn-2500") continue end
        # if (dir == "Sn-5000") continue end
        if (dir != "SNL-3-5000") && (restart == false) continue
        else restart = true end

        if maxrelax == 0.0
            outname = string("../OutfilesJO/Sn/",dir,"-", method, "-basic.out");
        else
            outname = string("../OutfilesJO/Sn/",dir,"-",method,"-",maxrelax,"-",lambda_reg,".out");
        end

        # scale the relaxation parameter as a function of the size of the instance
        maxrelax_scaled = maxrelax;
        if (dir == "Sn-100") maxrelax_scaled = sqrt(10.0) * maxrelax end
        if (dir == "Sn-500") maxrelax_scaled = sqrt(2.0) * maxrelax end
        if (dir == "Sn-5000") maxrelax_scaled = sqrt(0.2) * maxrelax end
        # outname = string("");

        run_directory(datasetpath, method, outname, maxrelax_scaled, lambda_reg, 0, norme, percent_closest);
    end
end

end # module
