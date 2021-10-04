using Printf

"""
Run the relaxed group transport with a range of relaxation parameters
bench_path: name of the benchmark directory
nbfiles: number of files considered per directory of the bench, 0 for all
norme : 1 or 2, norm used for distances in the space of covariates
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
                sol = ot_group(inst,0.2,maxrelax_group,norme,optimal);
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

