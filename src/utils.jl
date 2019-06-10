using JuMP
using Distances
using Printf
using DelimitedFiles
using DataFrames

@enum DataBase baseA baseB

# Definition and initialization of an Instance structure
struct Instance
    name::AbstractString
    nA::Int64
    nB::Int64
    Xobserv::Array{Float64,2}
    Yobserv::Array{Int64,1}
    Zobserv::Array{Int64,1}
    D::Array{Float64,2}
    Y::Array{Int64,1}
    Z::Array{Int64,1}
    indY::Dict{Int64,Array{Int64,1}}
    indZ::Dict{Int64,Array{Int64,1}}
    indXA:: Dict{Int64,Array{Int64}}; # indexes of subjects of A with given X value
    indXB:: Dict{Int64,Array{Int64}}; # indexes of subjects of B with given X value
    DA::Array{Float64,2}
    DB::Array{Float64,2}

  function Instance(data_file, norme::Int64)
      data = readdlm(data_file, ' ')

      # number of covariables
      nbcvar = size(data,2) - 3

      # recover the sets of individuals in base 1 and 2
      base = data[2:end,1]
      indA = findall(base .== 1)
      indB = findall(base .== 2)
      nA = length(indA)
      nB = length(indB)

      # recover the input data
      Xobserv = Array{Float64,2}(data[2:end, 4:end])
      Yobserv = Array{Float64,1}(data[2:end, 2])
      Zobserv = Array{Float64,1}(data[2:end, 3])

      # modify order so that base A comes first and then base B
      Xobserv = [Xobserv[indA,:];Xobserv[indB,:]];
      Yobserv = [Yobserv[indA];Yobserv[indB]];
      Zobserv = [Zobserv[indA];Zobserv[indB]];
      indA = 1:nA;
      indB = nA+1:nA+nB;

      # Modify Y and Z so that they go from 1 to the number of modalities
      Y = sort(unique(Yobserv[Yobserv .!= -1]));
      Z = sort(unique(Zobserv[Zobserv .!= -1]));
      for i in 1:length(Y)
          Yobserv[Yobserv .== Y[i]] .= i;
      end
      Y = [i for i in 1:length(Y)];
      for i in 1:length(Z)
          Zobserv[Zobserv .== Z[i]] .= i;
      end
      Z = [i for i in 1:length(Z)];

      # list the distinct modalities in A and B
      indY = Dict((m,findall(Yobserv[1:nA] .== m)) for m in Y)
      indZ = Dict((m,findall(Zobserv[nA+1:end] .== m)) for m in Z)

      # compute the distance between pairs of individuals in different bases
      # devectorize all the computations to go about twice faster
      # only compute norm 1 here
      a = transpose(Xobserv[indA,:])
      b = transpose(Xobserv[indB,:])
      if (norme == 1)
          D = pairwise(Cityblock(), a, b, dims=2)
          DA = pairwise(Cityblock(), a, a, dims=2)
          DB = pairwise(Cityblock(), b, b, dims=2)
      elseif (norme == 2)
          D = pairwise(Euclidean(), a, b, dims=2)
          DA = pairwise(Euclidean(), a, a, dims=2)
          DB = pairwise(Euclidean(), b, b, dims=2)
      elseif (norme == 0)
          D = pairwise(Hamming(), a, b, dims=2)
          DA = pairwise(Hamming(), a, a, dims=2)
          DB = pairwise(Hamming(), b, b, dims=2)
      end

      # Compute the indexes of individuals with same covariates
      A = 1:nA;
      B = 1:nB;
      nbX = 0;
      indXA = Dict{Int64,Array{Int64}}();
      indXB = Dict{Int64,Array{Int64}}();

      X1val = sort(unique(Xobserv[:,1]));
      X2val = sort(unique(Xobserv[:,2]));
      X3val = sort(unique(Xobserv[:,3]));

      Xval = convert(Matrix,unique(DataFrame(Xobserv)));

      # aggregate both bases
      for i in  1:size(Xval,1)
          nbX = nbX + 1;
          x = zeros(size(Xval,2),1);
          x[:,1] = [Xval[i,j] for j in 1:size(Xval,2)];
          if (norme == 1)
              distA = pairwise(Cityblock(), x, transpose(Xobserv[A,:]), dims=2);
              distB = pairwise(Cityblock(), x, transpose(Xobserv[B .+ nA,:]), dims=2);
          elseif (norme == 2)
              distA = pairwise(Euclidean(), x, transpose(Xobserv[A,:]), dims=2);
              distB = pairwise(Euclidean(), x, transpose(Xobserv[B .+ nA,:]), dims=2);
          elseif (norme == 0)
              distA = pairwise(Hamming(), x, transpose(Xobserv[A,:]), dims=2);
              distB = pairwise(Hamming(), x, transpose(Xobserv[B .+ nA,:]), dims=2);
          end
          indXA[nbX] = findall(distA[1,:] .< 0.1);
          indXB[nbX] = findall(distB[1,:] .< 0.1);
      end

      file_name = basename(data_file)
      new(file_name,nA, nB, Xobserv, Yobserv, Zobserv, D, Y, Z, indY, indZ, indXA, indXB, DA, DB)
  end
end


function aggregate_per_covar_mixed(inst::Instance, norme::Int64=1, aggregate_tol::Float64=0.5)

    A = 1:inst.nA;
    B = 1:inst.nB;

    # initialization of the structures
    nbX = 0;
    indXA = Dict{Int64,Array{Int64}}();
    indXB = Dict{Int64,Array{Int64}}();
    notaggA = [i for i in A];
    notaggB = [i for i in B];

    # aggregate until every individual in base A is aggregated
    while !isempty(notaggA)
        nbX += 1;
        ind = notaggA[1];
        isinset = inst.DA[ind,notaggA] .< aggregate_tol;
        indXA[nbX] = notaggA[isinset];
        deleteat!(notaggA, isinset);
        isinset = inst.D[ind,notaggB] .< aggregate_tol;
        indXB[nbX] = notaggB[isinset];
        deleteat!(notaggB, isinset);
    end

    # complete the aggregation with the individuals of base B that are not aggregated yet
    while !isempty(notaggB)
        nbX += 1;
        ind = notaggB[1];
        isinset = inst.DB[ind,notaggB] .< aggregate_tol;
        indXB[nbX] = notaggB[isinset];
        indXA[nbX] = [];
        deleteat!(notaggB, isinset);
    end

    return indXA, indXB;
end

###############################################################################
# Compute a bound on the average prediction error in each base
# The bound is computed as the expected prediction error assuming that the
# distribution of Z in base A (and that of Y in base B) is known, and the
# prediction done with the value that maximizes the probability
###############################################################################
function bound_prediction_error(inst::Instance, norme::Int64=1, aggregate_tol::Float64=0.5)

    # Local redefinitions of parameters of  the instance
    nA = inst.nA;
    nB = inst.nB;
    Y = copy(inst.Y);
    Z = copy(inst.Z);

    # compute the bound in base A
    boundpredZA = 0.0;
    for x  in 1:length(inst.indXA)
        for mA in Y
            indwithmA = inst.indXA[x][findall(inst.Yobserv[inst.indXA[x]] .== mA)];
            nindiv = length(indwithmA);
            if nindiv == 0
                continue
            end
            for mB in Z
                indwithmB = indwithmA[findall(inst.Zobserv[indwithmA] .== mB)];
                boundpredZA += length(indwithmB)/nindiv * (1-length(indwithmB)/nindiv) * nindiv/ nA;
            end
        end
    end
    # @printf("Bound on average prediction error in A : %.1f %%\n", 100.*boundpredZA);

    # compute the bound in base B
    boundpredYB = 0.0;
    for x  in 1:length(inst.indXB)
        for mB in Z
            indwithmB = inst.indXB[x][findall(inst.Zobserv[inst.indXB[x].+nA] .== mB)].+ nA;
            nindiv = length(indwithmB);
            if nindiv == 0
                continue
            end
            for mA in Y
                indwithmA = indwithmB[findall(inst.Yobserv[indwithmB] .== mA)];
                boundpredYB += length(indwithmA)/nindiv * (1-length(indwithmA)/nindiv) * nindiv/ nB;
            end
        end
    end
    # @printf("Bound on average prediction error in B : %.1f %%\n", 100.*boundpredYB);

    return boundpredZA, boundpredYB;
end

###############################################################################
# Return the empirical cardinality of the joint occurrences of (C=x,Y=mA,Z=mB)
# in both bases
###############################################################################
function empirical_distribution(inst::Instance, norme::Int64=0, aggregate_tol::Float64=0.5)

    # Local redefinitions of parameters of  the instance
    nA = inst.nA;
    nB = inst.nB;
    A = 1:inst.nA;
    B = 1:inst.nB;
    Y = copy(inst.Y);
    Z = copy(inst.Z);

    # aggregate the individuals per covariate
    nbX = length(inst.indXA);

    # count the cardinality of occurrence of each triplet (x,mA,mB) in both bases
    cardA_c_mA_mB = zeros(nbX, length(Y), length(Z));
    for x = 1:nbX
        for i in inst.indXA[x]
            cardA_c_mA_mB[x,inst.Yobserv[i],inst.Zobserv[i]] += 1;
        end
    end
    cardB_c_mA_mB = zeros(nbX, length(Y), length(Z));
    for x = 1:nbX
        for i in inst.indXB[x]
            cardB_c_mA_mB[x,inst.Yobserv[i+nA],inst.Zobserv[i+nA]] += 1;
        end
    end

    return cardA_c_mA_mB, cardB_c_mA_mB;
end


###############################################################################
# Display information about the distance between the modalities
###############################################################################

function disp_inst_info(inst::Instance)

    # local definitions
    nA = inst.nA
    nB = inst.nB
    A = 1:inst.nA
    B = 1:inst.nB
    Y = copy(inst.Y)
    Z = copy(inst.Z)
    indY = copy(inst.indY)
    indZ = copy(inst.indZ)


    println("\n#################################################################")
    println("INFORMATION ABOUT THE INSTANCE")
    println("#################################################################\n")

    # return indicators about the original density of the modalities
    print("Average distance between objects of base 1: ")
    @printf("%.2f\n", 10*sum([DA[i,j] for i in A, j in A])/(nA^2))
    print("Average distance between objects of base 2: ")
    @printf("%.2f\n", 10*sum([DB[i,j] for i in B, j in B])/(nB^2))
    print("Crossed average distance between objects of base 1 and 2: ")
    @printf("%.2f\n", 10*sum([inst.D[i,j] for i in A, j in B])/(nA*nB))

    # restrict the average distance to the 10% closest individuals
    println("\nAverage distance between objects per modality")
    println("Modalities of Y, individuals of A:")
    percent_closest = 0.1
    for y1 in Y
        for y2 in Y
            if y1 > y2
                continue
            end
            avg = avg_distance_closest(inst,baseA,baseA,baseA,y1,y2,1.0)
            @printf("\tModalities %d and %d : %.2f\n", y1, y2, 10*avg)
            avg = avg_distance_closest(inst,baseA,baseA,baseA,y1,y2,percent_closest)
            @printf("\t\trestricted to the %.1f %% closest: %.2f\n", 100*percent_closest, 10*avg)
        end
    end
    println("\nModalities of Z, individuals of B:")
    for z1 in Z
        for z2 in Z
            if z1 > z2
                continue
            end
            avg = avg_distance_closest(inst,baseB,baseB,baseB,z1,z2,1.0)
            @printf("\tModalities %d and %d : %.2f\n", z1, z2, 10*avg)
            avg = avg_distance_closest(inst,baseB,baseB,baseB,z1,z2,percent_closest)
            @printf("\t\trestricted to the %.1f %% closest: %.2f\n", 100*percent_closest, 10*avg)
        end
    end
    println("\nModalities of Y, crossed bases:")
    for y1 in Y
        for y2 in Y
            avg = avg_distance_closest(inst,baseA,baseB,baseA,y1,y2,1.0)
            @printf("\tModalities %d and %d : %.2f\n", y1, y2, 10*avg)
            avg = avg_distance_closest(inst,baseA,baseB,baseA,y1,y2,percent_closest)
            @printf("\t\trestricted to the %.1f %% closest: %.2f\n", 100*percent_closest, 10*avg)
        end
    end
    println("\nModalities of Z, crossed bases:")
    for z1 in Z
        for z2 in Z
            avg = avg_distance_closest(inst,baseA,baseB,baseB,z1,z2,1.0)
            @printf("\tModalities %d and %d : %.2f\n", z1, z2, 10*avg)
            avg = avg_distance_closest(inst,baseA,baseB,baseB,z1,z2,percent_closest)
            @printf("\t\trestricted to the %.1f%% closest: %.2f\n", 100.0*percent_closest, 10*avg)
        end
    end
end

###############################################################################
# Compute the average distance between individuals of base1 with modality m1
# for outcome and individuals of base2 with modality m2 for outcome
# Consider only the percent_closest individuals in the computation of the
# distance
###############################################################################

function avg_distance_closest(inst::Instance, base1::DataBase, base2::DataBase, outcome::DataBase, m1::Int, m2::Int, percent_closest::Float64)

    # indices of individuals of base A with given outcomes in base B (and reciprocally)
    indZinA = Dict((z,findall(inst.Zobserv[1:inst.nA] .== z)) for z in inst.Z)
    indYinB = Dict((m,findall(inst.Yobserv[inst.nA+1:end] .== y)) for y in inst.Y)

    ind1 = base1 == baseA ? (outcome == baseA ? inst.indY[m1] : indZinA[m1]) : (outcome == baseA ? indYinB[m1] : inst.indZ[m1])
    ind2 = base2 == baseA ? (outcome == baseA ? inst.indY[m2] : indZinA[m2]) : (outcome == baseA ? indYinB[m2] : inst.indZ[m2])

    # select the distance matrix depending on the base
    # println(base1)
    # println(base2)
    # println(outcome)
    # println(typeof(Dscaled))
    D = base1 == baseA ? ( base2==baseA ? inst.DA : inst.D) : ( base2==baseA ? inst.D : inst.DB)

    # swap the two sets of indices if base1=baseB and base2=baseA
    if (base1==baseB && base2==baseA)
        ind = ind1
        ind1 = ind2
        ind2 = ind1
    end

    # compute the average distance between the individuals in ind1 and the
    # percent_closest in ind2 and reciprocally
    avg = 0.0
    for i in ind1
        nbclose = round(Int,percent_closest*length(ind2))
        distance = sort([D[i,j] for j in ind2])
        avg += sum(distance[1:nbclose])/nbclose/length(ind1)/2
    end
    for j in ind2
        nbclose = round(Int,percent_closest*length(ind1))
        distance = sort([D[i,j] for i in ind1])
        avg += sum(distance[1:nbclose])/nbclose/length(ind2)/2
    end

    return avg
end



mutable struct Solution
    tsolve::Float64  # solution time
    jointYZA::Array{Float64,2} # joint distribution of Y and Z in A
    jointYZB::Array{Float64,2} # joint distribution of Y and Z in B
    estimatorZA::Array{Float64,3} # estimator of probability of Z for individuals in base A
    estimatorYB::Array{Float64,3} # estimator of probability of Y for individuals in base B
    errorpredZA::Float64
    errorpredYB::Float64
    errorpredavg::Float64
    errordistribZA::Float64
    errordistribYB::Float64
    errordistribavg::Float64

    Solution(t,jointYZA,jointYZB) = new(t,jointYZA,jointYZB);
    Solution(t,jointYZA,jointYZB,estimatorZA,estimatorYB) = new(t,jointYZA,jointYZB,estimatorZA,estimatorYB);
end


###############################################################################
# Compute prediction errors in a solution
###############################################################################
function compute_pred_error(inst::Instance, sol, proba_disp::Bool=true, mis_disp::Bool=false, full_disp::Bool=false)

    A = 1:inst.nA
    B = 1:inst.nB

    # display the transported and real modalities
    if full_disp
        println("Modalities of base 1 individuals:")
        for i in A
            println("Index: " , i, ", real value: ", inst.Zobserv[i], ", transported value: ", sol.predZA[i])
        end
        # display the transported and real modalities
        println("Modalities of base 2 individuals:")
        for j in B
            println("Index: " , j, ", real value: ", inst.Yobserv[inst.nA+j], ", transported value: ", sol.predYB[j])
        end
    end

    # Count the number of mistakes in the transport

    # Base 1
    nbmisA = 0
    misA = []
    for i in A
        if sol.predZA[i] != inst.Zobserv[i]
          nbmisA += 1
          push!(misA, i )
        end
    end

    # Base 2
    nbmisB = 0
    misB = []
    for j in B
        if sol.predYB[j] != inst.Yobserv[inst.nA+j]
          nbmisB += 1
          push!(misB, j)
        end
    end


    if proba_disp
        if nbmisA == 0
            println("No mistake in the transport of base A")
        else
            @printf("Probability of error in base A: %.1f %%\n", 100.0 * nbmisA / inst.nA);
            if mis_disp
                println("Indices with mistakes in base A:", misA)
            end
        end

        if nbmisB == 0
            println("No mistake in the transport of base B")
        else
            @printf("Probability of error in base B: %.1f %%\n", 100.0 * nbmisB/inst.nB);
            if mis_disp
                println("Indices with mistakes in base 2:", misB)
            end
        end
    end

    sol.errorpredZA, sol.errorpredYB = nbmisA/inst.nA, nbmisB/inst.nB;
    sol.errorpredavg = (inst.nA*sol.errorpredZA + inst.nB*sol.errorpredYB)/(inst.nA+inst.nB);

    return sol;
end

###############################################################################
# Compute errors in the conditional distributions of a solution
###############################################################################
function compute_distrib_error(inst::Instance, sol::Solution, empiricalZA, empiricalYB)

    nA = inst.nA ;  nB = inst.nB ; Y = copy(inst.Y); Z = copy(inst.Z);
    nbX = length(inst.indXA);


    sol.errordistribZA = sum(length(inst.indXA[x][findall(inst.Yobserv[inst.indXA[x]] .== y)])/nA * sum(max.(sol.estimatorZA[x,y,:] .- empiricalZA[x,y,:],0)) for x in 1:nbX, y in Y);
    sol.errordistribYB = sum(length(inst.indXB[x][findall(inst.Zobserv[inst.indXB[x].+nA] .== z)])/nB * sum(max.(sol.estimatorYB[x,:,z] .- empiricalYB[x,:,z],0)) for x in 1:nbX, z in Z);
    sol.errordistribavg = (nA * sol.errordistribZA + nB * sol.errordistribYB)/(nA+nB);

    return sol;
end

###############################################################################
# Compute the cost between pairs of outcomes as the average distance between
# covariations of individuals with these outcomes, but considering only the
# percent closest neighbors
###############################################################################

function average_distance_to_closest(inst::Instance, percent_closest::Float64)

    # Redefine A and B for the model
    A = 1:inst.nA
    B = 1:inst.nB
    Y = copy(inst.Y)
    Z = copy(inst.Z)
    indY = copy(inst.indY)
    indZ = copy(inst.indZ)

    # Compute average distances as described in the above
    Davg = zeros(length(Y),length(Z));
    DindivA  = zeros(inst.nA,length(Z));
    DindivB  = zeros(inst.nB,length(Y));
    for y in Y
        for i in indY[y]
            for z in Z
                nbclose = max(round(Int,percent_closest*length(indZ[z])),1);
                distance = sort([inst.D[i,j] for j in indZ[z]]);
                DindivA[i,z] =  sum(distance[1:nbclose])/nbclose;
                Davg[y,z] += sum(distance[1:nbclose])/nbclose/length(indY[y])/2.0;
            end
        end
    end
    for z in Z
        for j in indZ[z]
            for y in Y
                nbclose = max(round(Int,percent_closest*length(indY[y])),1);
                distance = sort([inst.D[i,j] for i in indY[y]]);
                DindivB[j,y] = sum(distance[1:nbclose])/nbclose;
                Davg[y,z] += sum(distance[1:nbclose])/nbclose/length(indZ[z])/2.0;
            end
        end
    end

    return Davg, DindivA, DindivB;
end
