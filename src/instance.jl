using DataFrames

export Instance

"""
    Instance(data_file, norme)

Definition and initialization of an Instance structure

- datafile : file name
- norme    : ( 1 : Cityblock, 2 : Euclidean, 3 : Hamming )
- indXA    : indexes of subjects of A with given X value
- indXB    : indexes of subjects of B with given X value
"""
struct Instance

    name    :: AbstractString
    nA      :: Int64
    nB      :: Int64
    Xobserv :: Array{Float64,2}
    Yobserv :: Array{Int64,1}
    Zobserv :: Array{Int64,1}
    D       :: Array{Float64,2}
    Y       :: Array{Int64,1}
    Z       :: Array{Int64,1}
    indY    :: Dict{Int64,Array{Int64,1}}
    indZ    :: Dict{Int64,Array{Int64,1}}
    indXA   :: Dict{Int64,Array{Int64}} 
    indXB   :: Dict{Int64,Array{Int64}} 
    DA      :: Array{Float64,2}
    DB      :: Array{Float64,2}

    function Instance(data_file :: String, norme::Int64)


      distance = Cityblock()

      if norme == 2
          distance = Euclidean()
      elseif norme == 3
          distance = Hamming()
      end

      data = readdlm(data_file, ' ')

      # number of covariables
      nbcvar = size(data,2) - 3

      # recover the sets of individuals in base 1 and 2
      base = data[2:end, 1]
      indA = findall(base .== 1)
      indB = findall(base .== 2)
      nA = length(indA)
      nB = length(indB)

      # recover the input data
      Xobserv = Array{Float64,2}(data[2:end, 4:end])
      Yobserv = Array{Float64,1}(data[2:end, 2])
      Zobserv = Array{Float64,1}(data[2:end, 3])

      # modify order so that base A comes first and then base B
      Xobserv = [Xobserv[indA, :]; Xobserv[indB, :]]
      Yobserv = [Yobserv[indA]; Yobserv[indB]]
      Zobserv = [Zobserv[indA]; Zobserv[indB]]
      indA = 1:nA
      indB = nA + 1:nA + nB

      # Modify Y and Z so that they go from 1 to the number of modalities
      Y = sort(unique(Yobserv[Yobserv .!= -1]))
      Z = sort(unique(Zobserv[Zobserv .!= -1]))
      for i in 1:length(Y)
          Yobserv[Yobserv .== Y[i]] .= i
      end
      Y = [i for i in 1:length(Y)]
      for i in 1:length(Z)
          Zobserv[Zobserv .== Z[i]] .= i
      end
      Z = [i for i in 1:length(Z)]

      # list the distinct modalities in A and B
      indY = Dict((m,findall(Yobserv[1:nA] .== m)) for m in Y)
      indZ = Dict((m,findall(Zobserv[nA + 1:end] .== m)) for m in Z)

      # compute the distance between pairs of individuals in different bases
      # devectorize all the computations to go about twice faster
      # only compute norm 1 here
      a = transpose(Xobserv[indA, :])
      b = transpose(Xobserv[indB, :])

      D = pairwise(distance, a, b, dims = 2)
      DA = pairwise(distance, a, a, dims = 2)
      DB = pairwise(distance, b, b, dims = 2)

      # Compute the indexes of individuals with same covariates
      A = 1:nA
      B = 1:nB
      nbX = 0
      indXA = Dict{Int64,Array{Int64}}()
      indXB = Dict{Int64,Array{Int64}}()

      X1val = sort(unique(Xobserv[:, 1]))
      X2val = sort(unique(Xobserv[:, 2]))
      X3val = sort(unique(Xobserv[:, 3]))

      Xval = convert(Matrix,unique(DataFrame(Xobserv)))

      # aggregate both bases
      x = zeros(size(Xval,2),1)
      for i in  1:size(Xval,1)
          nbX = nbX + 1
          fill!(x,0)
          x[:,1] .= Xval[i,:] 
          distA = pairwise(distance, x, transpose(Xobserv[A,:]), dims=2)
          distB = pairwise(distance, x, transpose(Xobserv[B .+ nA,:]), dims=2)
          indXA[nbX] = findall(distA[1,:] .< 0.1)
          indXB[nbX] = findall(distB[1,:] .< 0.1)
      end

      file_name = basename(data_file)
      new(file_name,nA, nB, Xobserv, Yobserv, Zobserv, 
          D, Y, Z, indY, indZ, indXA, indXB, DA, DB)
    end

    """
        Instance( df, distance)
    - df : dataframe with column names : ident,Y1,Y2,X1,X2,X3
    - distance : Cityblock(), Euclidean() or Hamming()
    """
    function Instance(df          :: DataFrame, 
                      covariables :: Vector{Symbol}, 
                      outcomes    :: Vector{Symbol}, 
                      distance    :: Distances.Metric )

      # number of covariables
      nbcvar = length(covariables)

      # recover the sets of individuals in base 1 and 2
      base = df.ident
      indA = findall(base .== 1)
      indB = findall(base .== 2)
      nA   = length(indA)
      nB   = length(indB)

      # recover the input data

      nobs = size(df,1)

      Xobserv = convert(Matrix, df[covariables])
      Yobserv = convert(Array, df.Y1)
      Zobserv = convert(Array, df.Y2)

      # modify order so that base A comes first and then base B
      Xobserv = [Xobserv[indA,:] ; Xobserv[indB,:]]
      Yobserv = [Yobserv[indA]   ; Yobserv[indB]  ]
      Zobserv = [Zobserv[indA]   ; Zobserv[indB]  ]

      indA = 1:nA
      indB = nA+1:nA+nB

      # Modify Y and Z so that they go from 1 to the number of modalities
      Y = sort(unique(Yobserv[Yobserv .!= -1]))
      Z = sort(unique(Zobserv[Zobserv .!= -1]))

      for i in 1:length(Y)
          Yobserv[Yobserv .== Y[i]] .= i
      end
      Y = eachindex(Y)
      for i in 1:length(Z)
          Zobserv[Zobserv .== Z[i]] .= i
      end
      Z = eachindex(Z)

      # list the distinct modalities in A and B
      indY = Dict((m,findall(Yobserv[1:nA]     .== m)) for m in Y)
      indZ = Dict((m,findall(Zobserv[nA+1:end] .== m)) for m in Z)

      # compute the distance between pairs of individuals in different bases
      # devectorize all the computations to go about twice faster
      # only compute norm 1 here
      a = transpose(Xobserv[indA,:])
      b = transpose(Xobserv[indB,:])

      D  = pairwise(distance, a, b, dims=2)
      DA = pairwise(distance, a, a, dims=2)
      DB = pairwise(distance, b, b, dims=2)

      # Compute the indexes of individuals with same covariates
      A     = 1:nA
      B     = 1:nB
      nbX   = 0
      indXA = Dict{Int64,Array{Int64}}()
      indXB = Dict{Int64,Array{Int64}}()

      X1val = sort(unique(Xobserv[:,1]))
      X2val = sort(unique(Xobserv[:,2]))
      X3val = sort(unique(Xobserv[:,3]))

      Xval  = convert(Matrix,unique(DataFrame(Xobserv)))

      # aggregate both bases
      for i in  1:size(Xval,1)
          nbX += 1
          x = zeros(size(Xval,2),1)
          x[:,1] = [Xval[i,j] for j in 1:size(Xval,2)]
          distA = pairwise(distance, x, transpose(Xobserv[A,:]), dims=2)
          distB = pairwise(distance, x, transpose(Xobserv[B .+ nA,:]), dims=2)
          indXA[nbX] = findall(distA[1,:] .< 0.1)
          indXB[nbX] = findall(distB[1,:] .< 0.1)
      end

      file_name = ""
      new(file_name,nA, nB, Xobserv, Yobserv, Zobserv, 
          D, Y, Z, indY, indZ, indXA, indXB, DA, DB)
  end

end
