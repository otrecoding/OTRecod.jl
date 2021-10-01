using MultivariateStats
using Distributions, PDMats


"""
    simulate


Simulate one dataset with three covariates described by their mean in each database (muA and muB) and the quantiles used for discretization (q1,q2,q3)
The dependency of outcomes on covariates is linear and given by the weights alpha1, alpha2 and by the R2 coefficient
The instance contains n individuals in each base
"""
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
    # XA = zeros(3,n)
    # XB = zeros(3,n)
    # XA[1,:] = [(UA[1,j] < t1[1] ? 1 : 2) for j in 1:n]
    # XA[2,:] = [(UA[2,j] < t2[1] ? 1 : (UA[2,j] < t2[2] ? 2 : 3)) for j in 1:n]
    # XA[3,:] = [(UA[3,j] < t3[1] ? 1 : (UA[3,j] < t3[2] ? 2 : (UA[3,j] < t3[3] ? 3 : 4))) for j in 1:n]
    # XB[1,:] = [(UB[1,j] < t1[1] ? 1 : 2) for j in 1:n]
    # XB[2,:] = [(UB[2,j] < t2[1] ? 1 : (UB[2,j] < t2[2] ? 2 : 3)) for j in 1:n]
    # XB[3,:] = [(UB[3,j] < t3[1] ? 1 : (UB[3,j] < t3[2] ? 2 : (UB[3,j] < t3[3] ? 3 : 4))) for j in 1:n]
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

    varerror = (1/R2 -1)*sum(Sigma);
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

