using Printf

# Write the dataset given by covariates and outcomes in input
function writedataset(outname,XA,YA,ZA,XB,YB,ZB)
    outfile = open(outname, "w");
    @printf(outfile,"id Y Z X1 X2 X3");
    n = length(YA);
    for i in 1:n
        @printf(outfile, "\n1 %d %d %d %d %d", YA[i], ZA[i], XA[1,i], XA[2,i], XA[3,i]);
    end
    for i in 1:n
        @printf(outfile, "\n2 %d %d %d %d %d", YB[i], ZB[i], XB[1,i], XB[2,i], XB[3,i]);
    end
    close(outfile);
end

# Write the dataset where some of the input covariates are considered as "latent" variables, meaning that outcomes depend on them, but that those covariates are not observed
# Instead of specifiying the latent variables, we provide the indices of those that are observed
function writedatasetlatent(outname,XA,YA,ZA,XB,YB,ZB,observed=[1,2,3])
    outfile = open(outname, "w");
    @printf(outfile,"id Y Z");
    for i in observed
        @printf(outfile," X%i", i);
    end
    n = length(YA);
    for i in 1:n
        @printf(outfile, "\n1 %d %d", YA[i], ZA[i]);
        for j in observed
            @printf(outfile, " %d", XA[j,i]);
        end
    end
    for i in 1:n
        @printf(outfile, "\n2 %d %d", YB[i], ZB[i]);
        for j in observed
            @printf(outfile, " %d", XB[j,i]);
        end
    end
    close(outfile);
end

