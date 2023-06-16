@enum DataBase baseA baseB

"""
$(SIGNATURES)

Compute the average distance between individuals of base1 with modality m1
for outcome and individuals of base2 with modality m2 for outcome

Consider only the percent_closest individuals in the computation of the
distance
"""
function avg_distance_closest(
    inst::Instance,
    base1::DataBase,
    base2::DataBase,
    outcome::DataBase,
    m1::Int,
    m2::Int,
    percent_closest::Float64,
)

    # indices of individuals of base A with given outcomes in base B (and reciprocally)
    indZinA = Dict((z, findall(inst.Zobserv[1:inst.nA] .== z)) for z in inst.Z)
    indYinB = Dict((m, findall(inst.Yobserv[inst.nA+1:end] .== y)) for y in inst.Y)

    ind1 =
        base1 == baseA ? (outcome == baseA ? inst.indY[m1] : indZinA[m1]) :
        (outcome == baseA ? indYinB[m1] : inst.indZ[m1])
    ind2 =
        base2 == baseA ? (outcome == baseA ? inst.indY[m2] : indZinA[m2]) :
        (outcome == baseA ? indYinB[m2] : inst.indZ[m2])

    # select the distance matrix depending on the base
    D =
        base1 == baseA ? (base2 == baseA ? inst.DA : inst.D) :
        (base2 == baseA ? inst.D : inst.DB)

    # swap the two sets of indices if base1=baseB and base2=baseA
    if (base1 == baseB && base2 == baseA)
        ind = ind1
        ind1 = ind2
        ind2 = ind1
    end

    # compute the average distance between the individuals in ind1 and the
    # percent_closest in ind2 and reciprocally
    avg = 0.0
    for i in ind1
        nbclose = round(Int, percent_closest * length(ind2))
        distance = sort([D[i, j] for j in ind2])
        avg += sum(distance[1:nbclose]) / nbclose / length(ind1) / 2
    end
    for j in ind2
        nbclose = round(Int, percent_closest * length(ind1))
        distance = sort([D[i, j] for i in ind1])
        avg += sum(distance[1:nbclose]) / nbclose / length(ind2) / 2
    end

    return avg
end


"""
$(SIGNATURES)

Compute the cost between pairs of outcomes as the average distance between
covariations of individuals with these outcomes, but considering only the
percent closest neighbors
"""
function average_distance_to_closest(inst::Instance, percent_closest::Float64)

    # Redefine A and B for the model
    A = 1:inst.nA
    B = 1:inst.nB
    Y = inst.Y
    Z = inst.Z
    indY = inst.indY
    indZ = inst.indZ

    # Compute average distances as described in the above
    Davg = zeros(Float64, (length(Y), length(Z)))
    DindivA = zeros(Float64, (inst.nA, length(Z)))
    DindivB = zeros(Float64, (inst.nB, length(Y)))

    for y in Y, i in indY[y], z in Z

        nbclose = max(round(Int, percent_closest * length(indZ[z])), 1)
        distance = sort([inst.D[i, j] for j in indZ[z]])
        DindivA[i, z] = sum(distance[1:nbclose]) / nbclose
        Davg[y, z] += sum(distance[1:nbclose]) / nbclose / length(indY[y]) / 2.0

    end

    for z in Z, j in indZ[z], y in Y

        nbclose = max(round(Int, percent_closest * length(indY[y])), 1)
        distance = sort([inst.D[i, j] for i in indY[y]])
        DindivB[j, y] = sum(distance[1:nbclose]) / nbclose
        Davg[y, z] += sum(distance[1:nbclose]) / nbclose / length(indZ[z]) / 2.0

    end

    Davg, DindivA, DindivB

end
