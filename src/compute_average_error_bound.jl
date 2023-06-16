"""
$(SIGNATURES)

Compute a lower bound on the best average prediction error that one can
obtain with a specific type of data sets
path: path of the directory containing the data set
"""
function compute_average_error_bound(path, norme::Int64 = 1)

    files = readdir(path)
    dir_name = split(path, "/")[end]
    errorboundA = 0.0
    errorboundB = 0.0
    nbfiles = 0

    for data_file in files
        # continue if not a data file
        if !(data_file[end-3:end] == ".txt")
            continue
        end

        # Reading the data file and preparing arrays
        inst = Instance(string(path, "/", data_file), norme)

        # compute the bound and update the cumulative value
        boundA, boundB = bound_prediction_error(inst, norme)
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

