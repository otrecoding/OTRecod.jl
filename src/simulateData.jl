using MultivariateStats
using Distributions, PDMats
using Printf



tabR2 = [0.01, 0.05, 0.1, 0.9]

# Simulate one dataset with three covariates described by their mean in each database (muA and muB) and the quantiles used for discretization (q1,q2,q3)
# The dependency of outcomes on covariates is linear and given by the weights alpha1, alpha2 and by the R2 coefficient
# The instance contains n individuals in each base
function simulate(R2     = 0.5,
                  muA    = [0.0, 0.0 ,0.0],
                  muB    = [1.0, 0.0, 0.0],
                  alphaA = [1.0  1.0  1.0],
                  alphaB = [1.0  1.0  1.0],
                  n = 1000,
                  q1 = [0.5],
                  q2 = [1.0/3.0, 2.0/3.0],
                  q3 = [0.25, 0.5, 0.75])

    Sigma = [[1.0,0.2,0.2] [0.2,1.0,0.2] [0.2,0.2,1.0]];
    PSigma = PDMat(Sigma);
    t1 = quantile.(Normal(muA[1], sqrt(Sigma[1,1])), q1 );
    t2 = quantile.(Normal(muA[1], sqrt(Sigma[1,1])), q2 );
    t3 = quantile.(Normal(muA[1], sqrt(Sigma[1,1])), q3 );
    UA = rand(MvNormal(muA, PSigma), n);
    UB = rand(MvNormal(muB, PSigma), n);

    XA = ones(3, n);
    XB = ones(3, n);
    for j in 1:n
        for i in 1:length(t1)
            if UA[1,j] >= t1[i] XA[1,j] = i+1; end
            if UB[1,j] >= t1[i] XB[1,j] = i+1; end
        end
        for i in 1:length(t2)
            if UA[2,j] >= t2[i] XA[2,j] = i+1; end
            if UB[2,j] >= t2[i] XB[2,j] = i+1; end
        end
        for i in 1:length(t3)
            if UA[3,j] >= t3[i] XA[3,j] = i+1; end
            if UB[3,j] >= t3[i] XB[3,j] = i+1; end
        end
    end

	dimU = 3;
    varerror = (1/R2 -1)*sum([alphaA[i]*alphaA[j]*Sigma[i,j] for i in 1:dimU, j in 1:dimU]);
    VA = transpose(alphaA * UA) .+ rand(Normal(0.0,sqrt(varerror)),n);
    VB = transpose(alphaB * UB) .+ rand(Normal(0.0,sqrt(varerror)),n);
    tY = quantile.(Normal((alphaA * muA)[1], sqrt(sum(alphaA * Sigma) + varerror)), [1.0/3.0, 2.0/3.0] );
    tZ = quantile.(Normal((alphaB * muB)[1], sqrt(sum(alphaB * Sigma) + varerror)), [0.25, 0.5, 0.75] );
    YA = [(VA[j] < tY[1] ? 1 : (VA[j] < tY[2] ? 2 : 3)) for j in 1:n];
    ZA = [(VA[j] < tZ[1] ? 1 : (VA[j] < tZ[2] ? 2 : (VA[j] < tZ[3] ? 3 : 4))) for j in 1:n];
    YB = [(VB[j] < tY[1] ? 1 : (VB[j] < tY[2] ? 2 : 3)) for j in 1:n];
    ZB = [(VB[j] < tZ[1] ? 1 : (VB[j] < tZ[2] ? 2 : (VB[j] < tZ[3] ? 3 : 4))) for j in 1:n];

    return XA, YA, ZA, XB, YB, ZB;
end

# Write the dataset given by covariates and outcomes in input
function writedataset(outname,XA,YA,ZA,XB,YB,ZB)
    outfile = open(outname, "w");
	n = length(YA);
	ncovar = size(XA)[1]
	@printf(outfile,"id Y Z");
	for i in 1:ncovar
		@printf(outfile," X%i", i);
	end
	for i in 1:n
        @printf(outfile, "\n1 %d %d", YA[i], ZA[i]);
        for j in 1:ncovar
            @printf(outfile, " %d", XA[j,i]);
        end
    end
    for i in 1:n
        @printf(outfile, "\n2 %d %d", YB[i], ZB[i]);
        for j in 1:ncovar
            @printf(outfile, " %d", XB[j,i]);
        end
    end
    close(outfile);
end


# Simulate default instance
for k = 1:100
    XA,YA,ZA,XB,YB,ZB = simulate()
    outname = "data/tab" * string(k) * ".txt"
    writedataset(outname,XA,YA,ZA,XB,YB,ZB)
end


# Simulate data with varying number of categories
mkpath( "../data/CAT-0")
mkpath( "../data/CAT-1")
mkpath( "../data/CAT-2")
mkpath( "../data/CAT-3")
for k = 1:100
    XA,YA,ZA,XB,YB,ZB = simulate(0.5,[0.0, 0.0 ,0.0],[1.0, 0.0, 0.0],[1.0 1.0 1.0],[1.0 1.0 1.0], 1000, [0.5], [1.0/3.0, 2.0/3.0],[0.5]);
    outname = "../data/CAT-0" * "/" * "tab" * string(k) * ".txt";
    writedataset(outname,XA,YA,ZA,XB,YB,ZB);
    XA,YA,ZA,XB,YB,ZB = simulate(0.5,[0.0, 0.0 ,0.0],[1.0, 0.0, 0.0],[1.0 1.0 1.0],[1.0 1.0 1.0], 1000, [0.5], [0.5],[0.5]);
    outname = "../data/CAT-1" * "/" * "tab" * string(k) * ".txt";
    writedataset(outname,XA,YA,ZA,XB,YB,ZB);
    XA,YA,ZA,XB,YB,ZB = simulate(0.5,[0.0, 0.0 ,0.0],[1.0, 0.0, 0.0],[1.0 1.0 1.0],[1.0 1.0 1.0], 1000, [1.0/3.0, 2.0/3.0], [1.0/3.0, 2.0/3.0], [1.0/3.0, 2.0/3.0]);
    outname = "../data/CAT-2" * "/" * "tab" * string(k) * ".txt";
    writedataset(outname,XA,YA,ZA,XB,YB,ZB);
    XA,YA,ZA,XB,YB,ZB = simulate(0.5,[0.0, 0.0 ,0.0],[1.0, 0.0, 0.0],[1.0 1.0 1.0],[1.0 1.0 1.0], 1000, [0.25, 0.5, 0.75], [0.25, 0.5, 0.75], [0.25, 0.5, 0.75]);
    outname = "../data/CAT-3" * "/" * "tab" * string(k) * ".txt";
    writedataset(outname,XA,YA,ZA,XB,YB,ZB);
end


# Simulate R2 instances
for k = 1:100
    for R2 in [0.01, 0.05, 0.1, 0.9]
        XA,YA,ZA,XB,YB,ZB = simulate(R2)
        path = "data/SR-" * string(R2)
        mkpath(path)
        outname = joinpath(path, "tab" * string(k) * ".txt")
        writedataset(outname,XA,YA,ZA,XB,YB,ZB)
    end
end

# Simulate alpha instances
mkpath( "data/Sa-1")
mkpath( "data/Sa-2")
mkpath( "data/Sa-3")

for k = 1:100
    XA,YA,ZA,XB,YB,ZB = simulate(0.5,[0.0, 0.0 ,0.0],[1.0, 0.0, 0.0],[1.0 1.0 1.0],[1.0 1.0 2.0])
    outname = "data/Sa-1" * "/" * "tab" * string(k) * ".txt"
    writedataset(outname,XA,YA,ZA,XB,YB,ZB)
    XA,YA,ZA,XB,YB,ZB = simulate(0.5,[0.0, 0.0 ,0.0],[1.0, 0.0, 0.0],[1.0 1.0 1.0],[1.0 1.5 2.0])
    outname = "data/Sa-2" * "/" * "tab" * string(k) * ".txt"
    writedataset(outname,XA,YA,ZA,XB,YB,ZB)
    XA,YA,ZA,XB,YB,ZB = simulate(0.5,[0.0, 0.0 ,0.0],[1.0, 0.0, 0.0],[1.0 1.0 1.0],[3.0 1.5 2.0])
    outname = "data/Sa-3" * "/" * "tab" * string(k) * ".txt"
    writedataset(outname,XA,YA,ZA,XB,YB,ZB)
end

# Simulate different distributions of the covariates
mkpath( "data/SX-1")
mkpath( "data/SX-2")
mkpath( "data/SX-3")
mkpath( "data/SX-4")
for k = 1:100
    XA,YA,ZA,XB,YB,ZB = simulate(0.5,[0.0, 0.0 ,0.0], [0.0, 0.0, 0.0])
    outname = joinpath("data/SX-1", "tab" * string(k) * ".txt")
    writedataset(outname,XA,YA,ZA,XB,YB,ZB)
    XA,YA,ZA,XB,YB,ZB = simulate(0.5,[0.0, 0.0 ,0.0], [0.5, 0.0, 0.0])
    outname = joinpath("data/SX-2", "tab" * string(k) * ".txt")
    writedataset(outname,XA,YA,ZA,XB,YB,ZB)
    XA,YA,ZA,XB,YB,ZB = simulate(0.5,[0.0, 0.0 ,0.0], [1.0, 1.0, 0.0])
    outname = joinpath("data/SX-3", "tab" * string(k) * ".txt")
    writedataset(outname,XA,YA,ZA,XB,YB,ZB)
    XA,YA,ZA,XB,YB,ZB = simulate(0.5,[0.0, 0.0 ,0.0], [1.0, 2.0, 0.0])
    outname = joinpath("data/SX-4", "tab" * string(k) * ".txt")
    writedataset(outname,XA,YA,ZA,XB,YB,ZB)
end

# Simulate with variations in the size of the databases
for k in 1:100
    for n in [250,2500]
        path = "data/Sn/Sn-" * string(n)
        mkpath(path)
        XA,YA,ZA,XB,YB,ZB = simulate(0.5,[0.0, 0.0 ,0.0],[1.0, 0.0, 0.0],[1.0 1.0 1.0],[1.0 1.0 1.0],n)
        outname = joinpath( path, "tab" * string(k) * ".txt")
        writedataset(outname,XA,YA,ZA,XB,YB,ZB)
    end
end

# Simulate heterogeneous groups
R2=0.5
muA = [0.0, 0.0 ,0.0]
muB = [1.0, 0.0, 0.0]
alphaA = [1.0 1.0 1.0]
alphaB = [1.0 1.0 1.0]
n=1000
path = "data/SHG"
mkpath(path)
for k in 1:100
    Sigma = [[1.0,0.2,0.2] [0.2,1.0,0.2] [0.2,0.2,1.0]]
    PSigma = PDMat(Sigma)
    t1 = quantile.(Normal(muA[1], sqrt(Sigma[1,1])), [0.5] )
    t2 = quantile.(Normal(muA[1], sqrt(Sigma[1,1])), [1.0/3.0, 2.0/3.0] )
    t3 = quantile.(Normal(muA[1], sqrt(Sigma[1,1])), [0.25, 0.5, 0.75] )
    UA = rand(MvNormal(muA, Sigma),n)
    UB = rand(MvNormal(muB, Sigma),n)
    XA = zeros(3,n)
    XB = zeros(3,n)
    XA[1,:] = [(UA[1,j] < t1[1] ? 1 : 2) for j in 1:n]
    XA[2,:] = [(UA[2,j] < t2[1] ? 1 : (UA[2,j] < t2[2] ? 2 : 3)) for j in 1:n]
    XA[3,:] = [(UA[3,j] < t3[1] ? 1 : (UA[3,j] < t3[2] ? 2 : (UA[3,j] < t3[3] ? 3 : 4))) for j in 1:n]
    XB[1,:] = [(UB[1,j] < t1[1] ? 1 : 2) for j in 1:n]
    XB[2,:] = [(UB[2,j] < t2[1] ? 1 : (UB[2,j] < t2[2] ? 2 : 3)) for j in 1:n]
    XB[3,:] = [(UB[3,j] < t3[1] ? 1 : (UB[3,j] < t3[2] ? 2 : (UB[3,j] < t3[3] ? 3 : 4))) for j in 1:n]

    varerror = (1/R2 -1)*sum(Sigma)
    VA = transpose(alphaA * UA) .+ rand(Normal(0.0,sqrt(varerror)),n)
    VB = transpose(alphaB * UB) .+ rand(Normal(0.0,sqrt(varerror)),n)
    tY = quantile.(Normal((alphaA * muA)[1], sqrt(sum(alphaA * Sigma) + varerror)), [1.0/3.0, 2.0/3.0] )
    tZ = quantile.(Normal((alphaB * muB)[1], sqrt(sum(alphaB * Sigma) + varerror)), [0.25, 0.5, 0.75] )
    YA = [(VA[j] < tY[1] ? 1 : (VA[j] < tY[2] ? 2 : 1)) for j in 1:n]
    ZA = [(VA[j] < tZ[1] ? 1 : (VA[j] < tZ[2] ? 2 : (VA[j] < tZ[3] ? 3 : 1))) for j in 1:n]
    YB = [(VB[j] < tY[1] ? 1 : (VB[j] < tY[2] ? 2 : 1)) for j in 1:n]
    ZB = [(VB[j] < tZ[1] ? 1 : (VB[j] < tZ[2] ? 2 : (VB[j] < tZ[3] ? 3 : 1))) for j in 1:n]

    outname = joinpath(path, "tab" * string(k) * ".txt")
    writedataset(outname,XA,YA,ZA,XB,YB,ZB)
end

# Simulate groups from two latent variables
path = "data/SHG-2"
mkpath(path)
for k in 1:100
    muA = [0.0, 0.0 ,0.0]
    muB = [1.0, 0.0, 0.0]
    Sigma = [[1.0,0.2,0.2] [0.2,1.0,0.2] [0.2,0.2,1.0]]
    PSigma = PDMat(Sigma)
    t1 = quantile.(Normal(muA[1], sqrt(Sigma[1,1])), [0.5] )
    t2 = quantile.(Normal(muA[1], sqrt(Sigma[1,1])), [1.0/3.0, 2.0/3.0] )
    t3 = quantile.(Normal(muA[1], sqrt(Sigma[1,1])), [0.25, 0.5, 0.75] )
    UA = rand(MvNormal(muA, Sigma),1000)
    UB = rand(MvNormal(muB, Sigma),1000)
    XA = zeros(3,1000)
    XB = zeros(3,1000)
    XA[1,:] = [(UA[1,j] < t1[1] ? 1 : 2) for j in 1:1000]
    XA[2,:] = [(UA[2,j] < t2[1] ? 1 : (UA[2,j] < t2[2] ? 2 : 3)) for j in 1:1000]
    XA[3,:] = [(UA[3,j] < t3[1] ? 1 : (UA[3,j] < t3[2] ? 2 : (UA[3,j] < t3[3] ? 3 : 4))) for j in 1:1000]
    XB[1,:] = [(UB[1,j] < t1[1] ? 1 : 2) for j in 1:1000]
    XB[2,:] = [(UB[2,j] < t2[1] ? 1 : (UB[2,j] < t2[2] ? 2 : 3)) for j in 1:1000]
    XB[3,:] = [(UB[3,j] < t3[1] ? 1 : (UB[3,j] < t3[2] ? 2 : (UB[3,j] < t3[3] ? 3 : 4))) for j in 1:1000]

    alpha1 = [2.0 1.0 0.0]
    alpha2 = [0.0 1.0 2.0]
    R2 = 0.5
    varerror = (1/R2 -1)*sum(Sigma)
    VA1 = transpose(alpha1 * UA) .+ rand(Normal(0.0,sqrt(varerror)),1000)
    VA2 = transpose(alpha2 * UA) .+ rand(Normal(0.0,sqrt(varerror)),1000)
    VB1 = transpose(alpha1 * UB) .+ rand(Normal(0.0,sqrt(varerror)),1000)
    VB2 = transpose(alpha2 * UB) .+ rand(Normal(0.0,sqrt(varerror)),1000)
    tY = quantile.(Chisq(2), [1.0/3.0, 2.0/3.0] )
    tZ = quantile.(Chisq(2), [0.25, 0.5, 0.75] )
    YA = [((VA1[j]/(sqrt(sum(Sigma)+varerror)))^2+(VA2[j]/(sqrt(sum(Sigma)+varerror)))^2 < tY[1] ? 1 : ((VA1[j]/(sqrt(sum(Sigma)+varerror)))^2+(VA2[j]/(sqrt(sum(Sigma)+varerror)))^2 < tY[2] ? 2 : 3)) for j in 1:1000]
    ZA = [((VA1[j]/(sqrt(sum(Sigma)+varerror)))^2+(VA2[j]/(sqrt(sum(Sigma)+varerror)))^2 < tZ[1] ? 1 : ((VA1[j]/(sqrt(sum(Sigma)+varerror)))^2+(VA2[j]/(sqrt(sum(Sigma)+varerror)))^2 < tZ[2] ? 2 : ((VA1[j]/(sqrt(sum(Sigma)+varerror)))^2+(VA2[j]/(sqrt(sum(Sigma)+varerror)))^2 < tZ[3] ? 3 : 4))) for j in 1:1000]
    YB = [( (VB1[j]/(sqrt(sum(Sigma)+varerror)))^2+ (VB2[j]/(sqrt(sum(Sigma)+varerror)))^2 < tY[1] ? 1 : ((VB1[j]/(sqrt(sum(Sigma)+varerror)))^2+ (VB2[j]/(sqrt(sum(Sigma)+varerror)))^2 < tY[2] ? 2 : 3)) for j in 1:1000]
    ZB = [((VB1[j]/(sqrt(sum(Sigma)+varerror)))^2+ (VB2[j]/(sqrt(sum(Sigma)+varerror)))^2 < tZ[1] ? 1 : ((VB1[j]/(sqrt(sum(Sigma)+varerror)))^2+ (VB2[j]/(sqrt(sum(Sigma)+varerror)))^2 < tZ[2] ? 2 : ((VB1[j] /(sqrt(sum(Sigma)+varerror)))^2+ (VB2[j]/(sqrt(sum(Sigma)+varerror)))^2 < tZ[3] ? 3 : 4))) for j in 1:1000]

    outname = joinpath( path, "tab" * string(k) * ".txt")
    outfile = open(outname, "w")
    @printf(outfile,"id Y Z X1 X2 X3")
    for i in 1:1000
        @printf(outfile, "\n1 %d %d %d %d %d", YA[i], ZA[i], XA[1,i], XA[2,i], XA[3,i])
    end
    for i in 1:1000
        @printf(outfile, "\n2 %d %d %d %d %d", YB[i], ZB[i], XB[1,i], XB[2,i], XB[3,i])
    end
    close(outfile)
end
