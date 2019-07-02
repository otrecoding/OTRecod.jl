"""
    Solution(t, jointYZA, jointYZB)

    Solution(t, jointYZA, jointYZB, estimatorZA, estimatorYB)

- tsolve       : solution time
- jointYZA     : joint distribution of Y and Z in A
- jointYZB     : joint distribution of Y and Z in B
- estimatorZA  : estimator of probability of Z for individuals in base A
- estimatorYB  : estimator of probability of Y for individuals in base B

"""
mutable struct Solution

    tsolve          :: Float64  
    jointYZA        :: Array{Float64,2} 
    jointYZB        :: Array{Float64,2}
    estimatorZA     :: Array{Float64,3}
    estimatorYB     :: Array{Float64,3}
    errorpredZA     :: Float64
    errorpredYB     :: Float64
    errorpredavg    :: Float64
    errordistribZA  :: Float64
    errordistribYB  :: Float64
    errordistribavg :: Float64

    Solution(t, jointYZA, jointYZB) = new(t,jointYZA,jointYZB)

    Solution(t, jointYZA, jointYZB, estimatorZA, estimatorYB) = new(t, 
             jointYZA, jointYZB, estimatorZA, estimatorYB)
end


