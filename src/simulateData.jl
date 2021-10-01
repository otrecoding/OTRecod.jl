using MultivariateStats
using Distributions, PDMats
using Printf



tabR2 = [0.01, 0.05, 0.1, 0.9]

# Simulate default instance
for k = 1:100
    XA,YA,ZA,XB,YB,ZB = simulate()
    outname = "data/tab" * string(k) * ".txt"
    writedataset(outname,XA,YA,ZA,XB,YB,ZB)
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


# Simulate instances with latent variables
mkpath( "../data/LA-0")
mkpath( "../data/LA-1")
mkpath( "../data/LA-2")
mkpath( "../data/LA-3")
for k = 1:100
    XA,YA,ZA,XB,YB,ZB = simulate(0.5,[0.0, 0.0 ,0.0],[1.0, 0.0, 0.0],[1.0 1.0 1.0],[1.0 1.0 1.0], 1000, [0.5], [1.0/3.0, 2.0/3.0],[0.5]);
    outname = "../data/LA-0" * "/" * "tab" * string(k) * ".txt";
    writedataset(outname,XA,YA,ZA,XB,YB,ZB);
    outname = "../data/LA-1" * "/" * "tab" * string(k) * ".txt";
    writedatasetlatent(outname,XA,YA,ZA,XB,YB,ZB,[2,3]);
    outname = "../data/LA-2" * "/" * "tab" * string(k) * ".txt";
    writedatasetlatent(outname,XA,YA,ZA,XB,YB,ZB,[1,3]);
    outname = "../data/LA-3" * "/" * "tab" * string(k) * ".txt";
    writedatasetlatent(outname,XA,YA,ZA,XB,YB,ZB,[1,2]);
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


# Plot all results
mkpath("OutfilesJO")
plot_scenario("OutfilesJO", "SX", "group-basic", 0.0, 0.7)
plot_scenario("OutfilesJO", "SX", "group-0.4-0.0", 0.0, 0.7)
plot_scenario("OutfilesJO", "SX", "joint-basic", 0.0, 0.7)
plot_scenario("OutfilesJO", "SX", "joint-0.4-0.1", 0.0, 0.7)

plot_scenario("OutfilesJO", "SR", "group-basic", 0.0, 0.7)
plot_scenario("OutfilesJO", "SR", "group-0.4-0.0", 0.0, 0.7)
plot_scenario("OutfilesJO", "SR", "joint-basic", 0.0, 0.7)
plot_scenario("OutfilesJO", "SR", "joint-0.4-0.1", 0.0, 0.7)

plot_scenario("OutfilesJO", "Sa", "group-basic", 0.0, 0.7)
plot_scenario("OutfilesJO", "Sa", "group-0.4-0.0", 0.0, 0.7)
plot_scenario("OutfilesJO", "Sa", "joint-basic", 0.0, 0.7)
plot_scenario("OutfilesJO", "Sa", "joint-0.4-0.1", 0.0, 0.7)

plot_scenario("OutfilesJO/Sn", "Sn", "group-basic", 0.0, 0.7)
plot_scenario("OutfilesJO/Sn", "Sn", "group-0.4-0.0", 0.0, 0.7)
plot_scenario("OutfilesJO/Sn", "Sn", "joint-basic", 0.0, 0.7)
plot_scenario("OutfilesJO/Sn", "Sn", "joint-0.4-0.1", 0.0, 0.7)

inst   = Instance("data/tab1.txt", 0)
x      = Array{Float64,2}(inst.Xobserv[1001:2000,:])
y      = Array{Float64,1}(inst.Yobserv[1001:2000])
mean_y = sum(y)/length(y)

sol    = MultivariateStats.llsq(x, y)
a, b   = sol[1:end-1], sol[end]
yp     = x * a .+ b

sstot  = sum((y .- mean_y).^2)
ssreg  = sum((yp .- mean_y).^2)
ssres  = sum((yp .- y).^2)

R2     = 1 - ssres/sstot

println("R2 = " , R2)
