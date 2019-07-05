"""
Run the transport of covariates and outcomes with a range of parameters for
relaxation and regularization
bench_path: name of the benchmark directory
nbfiles: number of files considered per directory of the bench, 0 for all
norme : 1 or 2, norm used for distances in the space of covariates
"""
function test_params_joint(bench_path, nbfiles::Int64=1, norme::Int64=1)

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
        if (dir == "Sn-500")  continue end
        if (dir == "Sn-50")   continue end
        if (dir == "Sn-100")  continue end
        if (dir == "S_NO")    continue end
        if (dir == "Spi-4")   continue end
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
                    sol = ot_joint(inst, maxrelax, lambda_reg, 0.2);
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

