using DelimitedFiles

export Instance

"""
$(TYPEDEF)

$(TYPEDSIGNATURES)

Definition and initialization of an Instance structure

- datafile : file name
- norme    : ( 1 : Cityblock, 2 : Euclidean, 3 : Hamming )
- indXA    : indexes of subjects of A with given X value
- indXB    : indexes of subjects of B with given X value
"""
struct Instance

    name::AbstractString
    nA::Int64
    nB::Int64
    Xobserv::Array{Float64,2}
    Yobserv::Array{Int64,1}
    Zobserv::Array{Int64,1}
    D::Array{Float64,2}
    Xval::Array{Float64,2}
    Y::Array{Int64,1}
    Z::Array{Int64,1}
    indY::Dict{Int64,Array{Int64,1}}
    indZ::Dict{Int64,Array{Int64,1}}
    indXA::Dict{Int64,Array{Int64}}
    indXB::Dict{Int64,Array{Int64}}
    DA::Array{Float64,2}
    DB::Array{Float64,2}

    function Instance(
        data_file::String,
        norme::Int64,
        observed::Array{Int64,1} = Array{Int64,1}(),
    )


        distance = Cityblock() # WeightedCityblock([1.0, 2.0, 3.0])

        if norme == 2
            distance = Euclidean()
        elseif norme == 0
            distance = Hamming() #WeightedHamming([2.0, 1.0])  #Hamming(); ##
        end

        data = readdlm(data_file, ' ')

        # number of covariables
        nbcvar = size(data, 2) - 3
        if isempty(observed)
            observed = Array(4:(4+nbcvar-1))
        else
            observed = observed .+ 3
        end

        # recover the sets of individuals in base 1 and 2
        base = data[2:end, 1]
        indA = findall(base .== 1)
        indB = findall(base .== 2)
        nA = length(indA)
        nB = length(indB)

        # recover the input data
        Xobserv = Array{Float64,2}(data[2:end, observed])
        Yobserv = Array{Float64,1}(data[2:end, 2])
        Zobserv = Array{Float64,1}(data[2:end, 3])

        # modify order so that base A comes first and then base B
        Xobserv = [Xobserv[indA, :]; Xobserv[indB, :]]
        Yobserv = [Yobserv[indA]; Yobserv[indB]]
        Zobserv = [Zobserv[indA]; Zobserv[indB]]
        indA = 1:nA
        indB = nA+1:nA+nB

        # Modify Y and Z so that they go from 1 to the number of modalities
        Y = sort(unique(Yobserv[Yobserv.!=-1]))
        Z = sort(unique(Zobserv[Zobserv.!=-1]))
        for i = 1:length(Y)
            Yobserv[Yobserv.==Y[i]] .= i
        end
        Y = [i for i = 1:length(Y)]
        for i = 1:length(Z)
            Zobserv[Zobserv.==Z[i]] .= i
        end
        Z = [i for i = 1:length(Z)]

        # list the distinct modalities in A and B
        indY = Dict((m, findall(Yobserv[1:nA] .== m)) for m in Y)
        indZ = Dict((m, findall(Zobserv[nA+1:end] .== m)) for m in Z)

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
        Xval = Matrix(unique(DataFrame(Xobserv, :auto)))

        # aggregate both bases
        x = zeros(size(Xval, 2), 1)
        for i = 1:size(Xval, 1)
            nbX = nbX + 1
            fill!(x, 0)
            x[:, 1] .= Xval[i, :]
            distA = pairwise(distance, x, transpose(Xobserv[A, :]), dims = 2)
            distB = pairwise(distance, x, transpose(Xobserv[B.+nA, :]), dims = 2)
            indXA[nbX] = findall(distA[1, :] .< 0.1)
            indXB[nbX] = findall(distB[1, :] .< 0.1)
        end

        file_name = basename(data_file)
        new(
            file_name,
            nA,
            nB,
            Xobserv,
            Yobserv,
            Zobserv,
            D,
            Xval,
            Y,
            Z,
            indY,
            indZ,
            indXA,
            indXB,
            DA,
            DB,
        )
    end

    """
        Instance( df, distance)

    - df : dataframe with column names : ident,Y1,Y2,X1,X2,X3
    - distance : Cityblock(), Euclidean() or Hamming()
    """
    function Instance(
        df::DataFrame,
        covariables::Vector{Symbol},
        outcomes::Vector{Symbol},
        distance::Distances.Metric,
    )

        # number of covariables
        nbcvar = length(covariables)

        # recover the sets of individuals in base 1 and 2
        base = df.ident
        indA = findall(base .== 1)
        indB = findall(base .== 2)
        nA = length(indA)
        nB = length(indB)

        # recover the input data

        nobs = size(df, 1)

        Xobserv = Matrix(df[!, covariables])
        Yobserv = convert(Array, df.Y1)
        Zobserv = convert(Array, df.Y2)

        # modify order so that base A comes first and then base B
        Xobserv = [Xobserv[indA, :]; Xobserv[indB, :]]
        Yobserv = [Yobserv[indA]; Yobserv[indB]]
        Zobserv = [Zobserv[indA]; Zobserv[indB]]

        indA = 1:nA
        indB = nA+1:nA+nB

        # Modify Y and Z so that they go from 1 to the number of modalities
        Y = sort(unique(Yobserv[Yobserv.!=-1]))
        Z = sort(unique(Zobserv[Zobserv.!=-1]))

        for i = 1:length(Y)
            Yobserv[Yobserv.==Y[i]] .= i
        end
        Y = eachindex(Y)
        for i = 1:length(Z)
            Zobserv[Zobserv.==Z[i]] .= i
        end
        Z = eachindex(Z)

        # list the distinct modalities in A and B
        indY = Dict((m, findall(Yobserv[1:nA] .== m)) for m in Y)
        indZ = Dict((m, findall(Zobserv[nA+1:end] .== m)) for m in Z)

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

        Xval = Matrix(unique(DataFrame(Xobserv, :auto)))

        # aggregate both bases
        for i = 1:size(Xval, 1)
            nbX += 1
            x = zeros(size(Xval, 2), 1)
            x[:, 1] = [Xval[i, j] for j = 1:size(Xval, 2)]
            distA = pairwise(distance, x, transpose(Xobserv[A, :]), dims = 2)
            distB = pairwise(distance, x, transpose(Xobserv[B.+nA, :]), dims = 2)
            indXA[nbX] = findall(distA[1, :] .< 0.1)
            indXB[nbX] = findall(distB[1, :] .< 0.1)
        end

        file_name = ""
        new(
            file_name,
            nA,
            nB,
            Xobserv,
            Yobserv,
            Zobserv,
            D,
            Xval,
            Y,
            Z,
            indY,
            indZ,
            indXA,
            indXB,
            DA,
            DB,
        )
    end

end
