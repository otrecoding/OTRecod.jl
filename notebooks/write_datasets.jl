using Printf

# Write the dataset given by covariates and outcomes in input
function writedataset(outname, XA, YA, ZA, XB, YB, ZB)
    outfile = open(outname, "w")
    n = length(YA)
    ncovar = size(XA)[1]
    @printf(outfile, "id Y Z")
    for i = 1:ncovar
        @printf(outfile, " X%i", i)
    end
    for i = 1:n
        @printf(outfile, "\n1 %d %d", YA[i], ZA[i])
        for j = 1:ncovar
            @printf(outfile, " %d", XA[j, i])
        end
    end
    for i = 1:n
        @printf(outfile, "\n2 %d %d", YB[i], ZB[i])
        for j = 1:ncovar
            @printf(outfile, " %d", XB[j, i])
        end
    end
    close(outfile)
end
