using OTRecod

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

function run_benchmark(
    path,
    method::Symbol,
    maxrelax::Float64 = 0.0,
    lambda_reg::Float64 = 0.0,
    norme::Int64 = 0,
    percent_closest::Float64 = 0.2,
    observed::Array{Int64,1} = Array{Int64,1}(),
)

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
        datasetpath = joinpath(path, dir)
        println(datasetpath)
        if !isdir(datasetpath)
            continue
        end
        # if (dir != "Sn-250") && (dir != "Sn-2500") continue end
        # if (dir == "Sn-5000") continue end
        # if (dir != "SNL-3-5000") && (restart == false) continue
        # else restart = true end

        if !isempty(observed)
            splitstr = split(dir, "-")
            dir = "SLA-" * splitstr[2]
        end
        if (maxrelax == 0.0) && (norme == 0)
            outname = string("../outfiles/", dir, "-", method, "-basic.out")
        else
            outname = string(
                "../outfiles/",
                dir,
                "-",
                method,
                "-",
                maxrelax,
                "-",
                lambda_reg,
                ".out",
            )
        end

        # scale the relaxation parameter as a function of the size of the instance
        maxrelax_scaled = maxrelax
        if (dir == "Sn-100")
            maxrelax_scaled = sqrt(10.0) * maxrelax
        end
        if (dir == "Sn-500")
            maxrelax_scaled = sqrt(2.0) * maxrelax
        end
        if (dir == "Sn-5000")
            maxrelax_scaled = sqrt(0.2) * maxrelax
        end


        run_directory(
            datasetpath,
            method,
            outname,
            maxrelax_scaled,
            lambda_reg,
            0,
            norme,
            percent_closest,
            observed,
        )
    end
end
