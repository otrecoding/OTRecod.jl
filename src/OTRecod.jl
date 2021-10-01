module OTRecod

using Statistics
using JuMP, Cbc, Clp
using Printf

include("utils.jl")
include("ot_group.jl")
include("ot_joint.jl")
include("plot_functions.jl")
include("simulate.jl")
include("write_datasets.jl")

export run_directory
export run_benchmark
export ot_group, ot_joint

"""
    compute_average_error_bound(path; norme=1)

Compute a lower bound on the best average prediction error that one can
obtain with a specific type of data sets
path: path of the directory containing the data set
"""
function compute_average_error_bound(path, norme::Int64=1)

    files = readdir(path)
    dir_name = split(path,"/")[end]
    errorboundA = 0.0
    errorboundB = 0.0
    nbfiles = 0

    for data_file in files
        # continue if not a data file
        if !(data_file[end - 3:end] == ".txt") continue end

        # Reading the data file and preparing arrays
        inst = Instance(string(path,"/",data_file), norme)

        # compute the bound and update the cumulative value
        boundA,boundB = bound_prediction_error(inst,norme)
        errorboundA += boundA
        errorboundB += boundB
        nbfiles += 1
    end

    errorboundA = errorboundA / nbfiles
    errorboundB = errorboundB / nbfiles

    @printf("Bound on average prediction error in A : %.1f %%\n", 100.0 * errorboundA)
    @printf("Bound on average prediction error in B : %.1f %%\n", 100.0 * errorboundB)

    errorboundA, errorboundB

end

"""
    empirical_estimator(path; observed = [])

Get an empirical estimator of the distribution of Z conditional to Y and X
on base A and reciprocally on base B obtain with a specific type of data sets

- `path`: path of the directory containing the data set
- `observed`: if nonempty, list of indices of the observed covariates; this allows to exclude some latent variables.

"""
function empirical_estimator(path,
                             observed   :: Array{Int64,1} = Array{Int64,1}())

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
    cardA_c_mA_mB = zeros(nbX,length(Y),length(Z))
    cardB_c_mA_mB = zeros(nbX,length(Y),length(Z))
    nbIndiv = 0
    for data_file in txt_files

        # read the data file and prepare arrays
        inst = Instance(joinpath(path,data_file), 0, observed)

        # update the cumulative count
        countA, countB  = empirical_distribution(inst, 0)
        cardA_c_mA_mB .+= countA
        cardB_c_mA_mB .+= countB
        nbIndiv += inst.nA
        if (nbIndiv >= 100000) break end

    end

    # Compute the empirical distribution of Z conditional to X and y in
    # base A and that of Y conditional to X and Z in base B
    # Z conditional to X and Y in base A
    empiricalZA = ones(nbX,length(Y),length(Z)) ./ length(Z)
    for x in 1:nbX, y in Y
        cardA_c_mA = sum(cardA_c_mA_mB[x, y, z] for z in Z)
        if cardA_c_mA > 0
            empiricalZA[x, y, :] .= cardA_c_mA_mB[x, y, :] ./ cardA_c_mA
        end
    end
    # Y conditional to X and Z in base B
    empiricalYB = ones(nbX,length(Y),length(Z)) ./ length(Y)
    for x in 1:nbX, z in Z
        cardB_c_mB = sum(cardB_c_mA_mB[x, y, z] for y in Y)
        if cardB_c_mB > 0
            empiricalYB[x, :, z] .= cardB_c_mA_mB[x, :, z] ./ cardB_c_mB
        end
    end

    empiricalZA, empiricalYB

end

"""
    run_directory(path, method; outname = "result.out",
                                maxrelax = 0.0,
                                lambda_reg = 0.0,
                                nbfiles = 0,
                                norme = 0,
                                percent_closest = 0.2,
                                observed = [])

 Run one given method on a given number of data files of a given directory
 The data files must be the only files with extension ".txt" in the directory

 - `path`   : name of the directory
 - `method` : `:group` or `:joint`
 - `maxrelax`: maximum percentage of deviation from expected probability masses
 - `lambda_reg`: coefficient measuring the importance of the regularization term
 - `nbfiles`: number of files considered, 0 if all the data files are tested
 - `norme`  : 0, 1 or 2, norm used for distances in the space of covariates
 - `percent_closest`: percent of closest neighbors taken in the computation of the costs (both distance and regularization related)
 - `observed`: if nonempty, list of indices of the observed covariates; this allows to exclude some latent variables.

"""

function run_directory(path            :: String,
                       method          :: Symbol,
                       outname         :: String  = "result.out",
                       maxrelax        :: Float64 = 0.0,
                       lambda_reg      :: Float64 = 0.0,
                       nbfiles         :: Int64   = 0,
                       norme           :: Int64   = 0,
                       percent_closest :: Float64 = 0.2,
                       observed        :: Array{Int64,1} = Array{Int64,1}())

    @info " RUN ONE METHOD ON THE FILES OF ONE DIRECTORY "
    @info " Method      : $method   "
    @info " Directory   : $path     "
    @info " Output file : $outname  "

    if (nbfiles > 0)
        println("\tTest only ", nbfiles, " files")
    else
        println("\tTest all the files of the directory")
    end

    println("\tRelaxation parameter: ", maxrelax)
    if (method == :joint)
       println("\tRegularization parameter: ", lambda_reg)
    end

    # initialize the output file
    outfile = open(outname, "w")
    @printf(outfile, "%-12s , %-6s , %-5s , %-5s , %-5s , %-6s , %-6s , %-8s , %-7s , %-7s , %-9s , %-6s\n","filename", "method", "norme", "relax", "regul", "ePredA", "ePredB", "ePredavg", "eDistrA", "eDistrB", "eDistravg", "cpu")
    close(outfile)

    # compute the true empirical conditional distributions for comparison with results
    @info "compute the empirical distributions of outcomes"
    empiricalZA,empiricalYB = empirical_estimator(path, observed)


    # compute a bound on the average prediction error in each base
    # println("... compute bounds on the average prediction errors\n")
    # errorboundA, errorboundB = compute_average_error_bound(path, norme)

    # solve the instances corresponding to each file
    files  = readdir(path)
    nbruns = 0
    for data_file in files
        # stop if the requested number of runs has been performed
        if ((nbfiles > 0) & (nbruns >= nbfiles)) break end
        # continue if not a data file
        if !(data_file[end-3:end]==".txt") continue end

        # Reading the data file and preparing arrays
        inst  = Instance(string(path,"/",data_file), norme, observed)
        nA    = inst.nA
        nB    = inst.nB
        Y     = inst.Y
        Z     = inst.Z
        indXA = inst.indXA
        indXB = inst.indXB
        nbX   = length(indXA)

        @info " File : $(joinpath(path,data_file)) "
        if method == :group
            indiv_method = maxrelax > 0.0 ? :optimal : :sequential
            sol = ot_group(inst,percent_closest,maxrelax,indiv_method)
            #PN lambda_reg = 0.0
        elseif method == :joint
            sol = ot_joint(inst, maxrelax, lambda_reg, percent_closest)
        end

        sol = compute_pred_error(inst, sol, false)
        sol = compute_distrib_error(inst, sol, empiricalZA, empiricalYB)
        @show sol
        # if (size(inst.Xval)[2] >= 3)
        #     sol = compute_distrib_error_3covar(sol, inst, empiricalZA, empiricalYB)
        #     @show sol
        # end

        outfile = open(outname, "a")
        @printf(outfile, "%-12s , %-6s , %-5d , %-5.2f , %-5.2f , %-6.3f , %-6.3f , %-8.3f , %-7.3f , %-7.3f , %-9.3f , %-6.2f\n", data_file, method, norme, maxrelax, lambda_reg, sol.errorpredZA, sol.errorpredYB, sol.errorpredavg, sol.errordistribZA, sol.errordistribYB, sol.errordistribavg, sol.tsolve)
        close(outfile)
        nbruns += 1
    end
end



"""
    run_benchmark(path, method;
                        maxrelax = 0.0,
                        lambda_reg = 0.0,
                        norme = 0,
                        percent_closest = 0.2,
                        observed = [])

 Run one method on the complete benchmark
 - `path`   : name of the benchmark directory
 - `method` : `:group` or `:joint`
 - `maxrelax`: maximum percentage of deviation from expected probability masses
 - `lambda_reg`: coefficient measuring the importance of the regularization term
 - `norme`  : 0, 1 or 2, norm used for distances in the space of covariates
 - `percent_closest`: percent of closest neighbors taken in the computation of the costs (both distance and regularization related)
 - `observed`: if nonempty, list of indices of the observed covariates; this allows to exclude some latent variables.

"""

function run_benchmark(path,
                       method           :: Symbol,
                       maxrelax         :: Float64 = 0.0,
                       lambda_reg       :: Float64 = 0.0,
                       norme            :: Int64 = 0,
                       percent_closest  :: Float64 = 0.2,
                       observed         :: Array{Int64,1} = Array{Int64,1}())

    println("RUN ONE METHOD ON THE COMPLETE BENCHMARK ")
    println("\tMethod: ", method)
    println("\tRelaxation parameter: ", maxrelax)
    if (method == :joint)
       println("\tRegularization parameter: ", lambda_reg)
    end

    dirlist = readdir(path)
    println(dirlist)
    restart = true
    for dir in dirlist
        datasetpath = joinpath(path,dir)
        println(datasetpath)
        if !isdir(datasetpath) continue end
        # if (dir != "Sn-250") && (dir != "Sn-2500") continue end
        # if (dir == "Sn-5000") continue end
        # if (dir != "SNL-3-5000") && (restart == false) continue
        # else restart = true end

        if !isempty(observed)
            splitstr = split(dir, "-");
            dir = "SLA-" * splitstr[2];
        end
        if (maxrelax == 0.0) && (norme == 0)
            outname = string("../outfiles/",dir,"-", method, "-basic.out");
        else
            outname = string("../outfiles/",dir,"-",method,"-",maxrelax,"-",lambda_reg, ".out");
        end

        # scale the relaxation parameter as a function of the size of the instance
        maxrelax_scaled = maxrelax
        if (dir == "Sn-100")  maxrelax_scaled = sqrt(10.0) * maxrelax end
        if (dir == "Sn-500")  maxrelax_scaled = sqrt(2.0)  * maxrelax end
        if (dir == "Sn-5000") maxrelax_scaled = sqrt(0.2)  * maxrelax end


        run_directory(datasetpath, method, outname, maxrelax_scaled,
                      lambda_reg, 0, norme, percent_closest, observed)
    end
end



"""
    test_params_group(bench_path, nbfiles::Int64=1, norme::Int64=0)

Run the relaxed group transport with a range of relaxation parameters
- `bench_path`: name of the benchmark directory
- `nbfiles`: number of files considered per directory of the bench, 0 for all
-  `norme` : 0, 1 or 2, norm used for distances in the space of covariates
"""
function test_params_group(bench_path, nbfiles::Int64=1, norme::Int64=1)

    println("\n#################################################################")
    println("TEST THE RELAXATION PARAMETER ON THE WHOLE BENCHMARK ")
    println("\n#################################################################\n")

    res_file_name = "params_group.out"
    println("Results written in ", res_file_name)
    println("Open maxrelax_group_template.out for the descriptions of the entries")
    resfile = open(res_file_name, "w");
    @printf(resfile, "%-12s , %-4s , %-5s , %-6s , %-6s\n","dataname", "norm", "relax", "error", "cpu");
    close(resfile);


    dir_names = readdir(bench_path)
    nbruns = 0
    for dir in dir_names
        if (!isdir(string(bench_path,"/",dir))) continue end

        # Test only on the instances with 1000 individuals per base
        if (dir == "Sn-5000") continue end
        if (dir == "Sn-500") continue end
        if (dir == "Sn-50") continue end
        if (dir == "Sn-100") continue end
        if (dir == "S_NO") continue end

        println("\nDIRECTORY: ", dir);
        empiricalZA,empiricalYB = empirical_estimator(string(bench_path,"/",dir), norme);
        file_names = readdir(string(bench_path,"/",dir));
        for maxrelax_group in 0.0:0.1:1.0
            println("\n#################################################################")
            println("MAXRELAX =  ", maxrelax_group)
            println("\n#################################################################\n")

            nbruns = 0 ;
            error_avg = 0.;
            tsolve = 0.;
            for data_file in file_names
                # stop if the requested number of runs has been performed
                if ((nbfiles > 0) & (nbruns >= nbfiles)) break end
                # continue if not a data file
                if !(data_file[end-3:end]==".txt") continue end

                # Reading the data file and preparing arrays
                inst = Instance(string(bench_path,"/",dir,"/",data_file), norme)

                # Run the group transport
                sol = OT_group(inst,0.2,maxrelax_group,norme,optimal);
                sol = compute_pred_error(inst, sol, false);
                sol = compute_distrib_error(inst, sol, empiricalZA, empiricalYB);
                error_avg += sol.errordistribavg/nbfiles;
                tsolve += sol.tsolve;
                nbruns += 1;
            end

            resfile = open(res_file_name, "a");
            @printf(resfile, "%-12s , %-4d , %-5.2f , %-6.3f , %-6.2f\n", dir, norme, maxrelax_group, error_avg, tsolve/nbruns);
            close(resfile);
        end
    end
end


"""
    test_params_joint(bench_path, nbfiles::Int64=1, norme::Int64=0)

Run the transport of covariates and outcomes with a range of parameters for relaxation and regularization
- `bench_path`: name of the benchmark directory
- `nbfiles`: number of files considered per directory of the bench, 0 for all
-  `norme` : 0, 1 or 2, norm used for distances in the space of covariates
"""
function test_params_joint(bench_path, nbfiles::Int64=1, norme::Int64=0)

    println("\n#################################################################")
    println("TEST THE PARAMETERS OF THE JOINT TRANSPORT ON THE WHOLE BENCHMARK ")
    println("\n#################################################################\n")

    res_file_name = "params_test_joint.out"
    println("Results written in ", res_file_name)
    # resfile = open(res_file_name, "w");
    # @printf(resfile, "%-12s , %-4s , %-5s , %-5s , %-6s , %-6s\n","dataname", "norm", "relax", "regul", "error", "cpu");
    # close(resfile);


    dir_names = readdir(bench_path)
    nbruns = 0
    restart = true
    for dir in dir_names

        if (!isdir(string(bench_path,"/",dir))) continue end
        println(dir)
        # Test only on the instances with 1000 individuals per base
        if (dir == "Sn-5000") continue end
        if (dir == "Sn-500") continue end
        if (dir == "Sn-50") continue end
        if (dir == "Sn-100") continue end
        if (dir == "S_NO") continue end
        if (dir == "Spi-4") continue end
        if (dir != "Spi-4") && (restart == false) continue
        else restart = true end

        println("\nDIRECTORY: ", dir);
        empiricalZA,empiricalYB = empirical_estimator(string(bench_path,"/",dir), norme);
        file_names = readdir(string(bench_path,"/",dir));
        reg_range = [0.0 0.01 0.05 0.1 0.5 1.0 5.0 10.0];
        for maxrelax in 0.0:0.2:1.0
            for lambda_reg in reg_range
                println("\n#################################################################")
                println("MAXRELAX =  ", maxrelax)
                println("LAMBDA_REG = ", lambda_reg)
                println("\n#################################################################\n")

                nbruns = 0
                error_avg = 0.;
                tsolve = 0.;
                for data_file in file_names
                    # stop if the requested number of runs has been performed
                    if ((nbfiles > 0) & (nbruns >= nbfiles)) break end
                    # continue if not a data file
                    if !(data_file[end-3:end]==".txt") continue end

                    # Reading the data file and preparing arrays
                    inst = Instance(string(bench_path,"/",dir,"/",data_file), norme)

                    # Run the group transport
                    sol = OT_joint(inst, maxrelax, lambda_reg, 0.2);
                    sol = compute_pred_error(inst, sol, false);
                    sol = compute_distrib_error(inst, sol, empiricalZA, empiricalYB);
                    error_avg += sol.errordistribavg/nbfiles;
                    tsolve += sol.tsolve;
                    nbruns += 1;
                end

                resfile = open(res_file_name, "a");
                @printf(resfile, "%-12s , %-4d , %-5.2f , %-5.2f , %-6.3f , %-6.2f\n", dir, norme, maxrelax, lambda_reg, error_avg, tsolve/nbruns);
                close(resfile);
            end
        end
    end
end


"""
    test_ncds(file_name)

Run every method on an input file generated from the real NCDS dataset
"""
function test_ncds(file_name)
    inst= Instance(file_name,0)
    solgroup=OT_group(inst)
    solgrouprelax=OT_group(inst,0.2, 0.199)
    soljoint=OT_joint(inst, 0.0, 0.0, 0.2)
    soljointrelax=OT_joint(inst, 0.199, 0.1, 0.2)
    effjoint = round.(Int,inst.nB*soljoint.jointYZB+inst.nA*soljoint.jointYZA);
    effjointrelax = round.(Int,inst.nB*soljointrelax.jointYZB+inst.nA*soljointrelax.jointYZA);
    effgroup = round.(Int,inst.nB*solgroup.jointYZB+inst.nA*solgroup.jointYZA);
    effgrouprelax = round.(Int,inst.nB*solgrouprelax.jointYZB+inst.nA*solgrouprelax.jointYZA);
    tab = [sum((inst.Yobserv.== y) .& (inst.Zobserv.==z)) for y in inst.Y, z in inst.Z];
    println("Error with group = ", 1/2.0*sum(abs.(tab.-effgroup))/sum(tab));
    println("Error with relaxed group = ",1/2.0*sum(abs.(tab.-effgrouprelax))/sum(tab));
    println("Error with  joint = ",1/2.0*sum(abs.(tab.-effjoint))/sum(tab));
    println("Error with relaxed joint = ",1/2.0*sum(abs.(tab.-effjointrelax))/sum(tab));
    return tab,effgroup,effgrouprelax,effjoint,effjointrelax
end

end # module
