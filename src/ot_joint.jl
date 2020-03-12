export ot_joint

"""
Model where we directly compute the distribution of the outcomes for each
individual or for sets of indviduals that similar values of covariates

- aggregate_tol: quantify how much individuals' covariates must be close for aggregation
- reg_norm: norm1, norm2 or entropy depending on the type of regularization
- percent_closest: percent of closest neighbors taken into consideration in regularization
- lambda_reg: coefficient measuing the importance of the regularization term
- full_disp: if true, write the transported value of each individual; otherwise, juste write the number of missed transports
- solver_disp: if false, do not display the outputs of the solver
"""
function ot_joint(inst            :: Instance,
                  maxrelax        :: Float64 = 0.0,
                  lambda_reg      :: Float64 = 0.0,
                  percent_closest :: Float64 = 0.2,
                  norme           :: Metric  = Cityblock(),
                  aggregate_tol   :: Float64 = 0.5,
                  full_disp       :: Bool    = false,
                  solver_disp     :: Bool    = false)

    @info " AGGREGATE INDIVIDUALS WRT COVARIATES               "
    @info " Reg. weight           = $(lambda_reg)              "
    @info " Percent closest       = $(100.0 * percent_closest) % "
    @info " Aggregation tolerance = $(aggregate_tol)           "

    tstart = time()

    # Local redefinitions of parameters of  the instance
    nA = inst.nA
    nB = inst.nB
    A = 1:nA
    B = 1:nB
    Y = inst.Y
    Z = inst.Z
    indY = inst.indY
    indZ = inst.indZ
    Xobserv = inst.Xobserv
    Yobserv = inst.Yobserv
    Zobserv = inst.Zobserv

    # Create a model for the optimal transport of individuals
    modelA = Model(with_optimizer(Clp.Optimizer,LogLevel = 0))
    modelB = Model(with_optimizer(Clp.Optimizer,LogLevel = 0))

    # Compute data for aggregation of the individuals
    # println("... aggregating individuals")
    indXA = inst.indXA
    indXB = inst.indXB
    nbX = length(indXA)

    # compute the neighbors of the covariates for regularization
    Xvalues = convert(Matrix,unique(DataFrame(Xobserv)))
    dist_X = pairwise(norme, transpose(Xvalues),transpose(Xvalues), dims = 2)
    voisins_X = dist_X .<= 1

    # println("... computing costs")
    C = average_distance_to_closest(inst, percent_closest)[1]

    # Compute the estimators that appear in the model

    estim_XA = Dict([(x, length(indXA[x]) / nA) for x in 1:nbX])
    estim_XB = Dict([(x, length(indXB[x]) / nB) for x in 1:nbX])
    estim_XA_YA = Dict([((x,y), length(indXA[x][findall(Yobserv[indXA[x]] .== y)]) / nA) for x in 1:nbX, y in Y])
    estim_XB_ZB = Dict([((x,z), length(indXB[x][findall(Zobserv[indXB[x].+nA] .== z)]) / nB) for x in 1:nbX, z in Z])


    # Basic part of the model

    # Variables
    # - gammaA[x,y,z]: joint probability of X=x, Y=y and Z=z in base A
    @variable(modelA, gammaA[x in 1:nbX, y in Y, z in Z] >= 0, base_name="gammaA")

    # - gammaB[x,y,z]: joint probability of X=x, Y=y and Z=z in base B
    @variable(modelB, gammaB[x in 1:nbX, y in Y, z in Z] >= 0, base_name="gammaB")

    @variable(modelA, errorA_XY[x in 1:nbX, y in Y], base_name="errorA_XY")
    @variable(modelA, abserrorA_XY[x in 1:nbX, y in Y] >= 0, base_name="abserrorA_XY")
    @variable(modelA, errorA_XZ[x in 1:nbX, z in Z], base_name="errorA_XZ")
    @variable(modelA, abserrorA_XZ[x in 1:nbX, z in Z] >= 0, base_name="abserrorA_XZ")

    @variable(modelB, errorB_XY[x in 1:nbX, y in Y], base_name="errorB_XY")
    @variable(modelB, abserrorB_XY[x in 1:nbX, y in Y] >= 0, base_name="abserrorB_XY")
    @variable(modelB, errorB_XZ[x in 1:nbX, z in Z], base_name="errorB_XZ")
    @variable(modelB, abserrorB_XZ[x in 1:nbX, z in Z] >= 0, base_name="abserrorB_XZ")

    # Constraints
    # - assign sufficient probability to each class of covariates with the same outcome
    @constraint(modelA, ctYandXinA[x in 1:nbX, y in Y], sum(gammaA[x,y,z] for z in Z) == estim_XA_YA[x,y] + errorA_XY[x,y])
    @constraint(modelB, ctZandXinB[x in 1:nbX, z in Z], sum(gammaB[x,y,z] for y in Y) == estim_XB_ZB[x,z] + errorB_XZ[x,z])

    # - we impose that the probability of Y conditional to X is the same in the two databases
    # - the consequence is that the probability of Y and Z conditional to Y is also the same in the two bases
    @constraint(modelA, ctZandXinA[x in 1:nbX, z in Z], estim_XB[x]*sum(gammaA[x,y,z] for y in Y) == estim_XB_ZB[x,z] * estim_XA[x] + estim_XB[x]*errorA_XZ[x,z])

    @constraint(modelB, ctYandXinB[x in 1:nbX, y in Y], estim_XA[x]*sum(gammaB[x,y,z] for z in Z) == estim_XA_YA[x,y] * estim_XB[x] + estim_XA[x] * errorB_XY[x,y])

    # - recover the norm 1 of the error
    @constraint(modelA,[x in 1:nbX, y in Y], errorA_XY[x,y] <= abserrorA_XY[x,y])
    @constraint(modelA,[x in 1:nbX, y in Y], -errorA_XY[x,y] <= abserrorA_XY[x,y])
    @constraint(modelA, sum( abserrorA_XY[x,y] for x in 1:nbX, y in Y) <= maxrelax/2.0)
    @constraint(modelA, sum(errorA_XY[x,y] for x in 1:nbX, y in Y) == 0.0)
    @constraint(modelA,[x in 1:nbX, z in Z], errorA_XZ[x,z] <= abserrorA_XZ[x,z])
    @constraint(modelA,[x in 1:nbX, z in Z], -errorA_XZ[x,z] <= abserrorA_XZ[x,z])
    @constraint(modelA, sum( abserrorA_XZ[x,z] for x in 1:nbX, z in Z) <= maxrelax/2.0)
    @constraint(modelA, sum(errorA_XZ[x,z] for x in 1:nbX, z in Z) == 0.0)

    @constraint(modelB,[x in 1:nbX, y in Y], errorB_XY[x,y] <= abserrorB_XY[x,y])
    @constraint(modelB,[x in 1:nbX, y in Y], -errorB_XY[x,y] <= abserrorB_XY[x,y])
    @constraint(modelB, sum( abserrorB_XY[x,y] for x in 1:nbX, y in Y) <= maxrelax/2.0)
    @constraint(modelB, sum(errorB_XY[x,y] for x in 1:nbX, y in Y) == 0.0)
    @constraint(modelB,[x in 1:nbX, z in Z], errorB_XZ[x,z] <= abserrorB_XZ[x,z])
    @constraint(modelB,[x in 1:nbX, z in Z], -errorB_XZ[x,z] <= abserrorB_XZ[x,z])
    @constraint(modelB, sum( abserrorB_XZ[x,z] for x in 1:nbX, z in Z) <= maxrelax/2.0)
    @constraint(modelB, sum(errorB_XZ[x,z] for x in 1:nbX, z in Z) == 0.0)

    # - regularization
    @variable(modelA, reg_absA[x1 in 1:nbX, x2 in findall(voisins_X[x1,:]), y in Y, z in Z] >= 0)
    @constraint(modelA, [x1 in 1:nbX, x2 in findall(voisins_X[x1,:]), y in Y, z in Z], reg_absA[x1,x2,y,z] >= gammaA[x1,y,z]/(max(1,length(indXA[x1]))/nA)-gammaA[x2,y,z]/(max(1,length(indXA[x2]))/nA))
    @constraint(modelA, [x1 in 1:nbX, x2 in findall(voisins_X[x1,:]), y in Y, z in Z], reg_absA[x1,x2,y,z] >= gammaA[x2,y,z]/(max(1,length(indXA[x2]))/nA)-gammaA[x1,y,z]/(max(1,length(indXA[x1]))/nA))
    @expression(modelA, regterm, sum(1/length(voisins_X[x1,:]) *reg_absA[x1,x2,y,z] for x1 in 1:nbX, x2 in findall(voisins_X[x1,:]), y in Y, z in Z))

    @variable(modelB, reg_absB[x1 in 1:nbX, x2 in findall(voisins_X[x1,:]), y in Y, z in Z] >= 0)
    @constraint(modelB, [x1 in 1:nbX, x2 in findall(voisins_X[x1,:]), y in Y, z in Z], reg_absB[x1,x2,y,z] >= gammaB[x1,y,z]/(max(1,length(indXB[x1]))/nB)-gammaB[x2,y,z]/(max(1,length(indXB[x2]))/nB))
    @constraint(modelB, [x1 in 1:nbX, x2 in findall(voisins_X[x1,:]), y in Y, z in Z], reg_absB[x1,x2,y,z] >= gammaB[x2,y,z]/(max(1,length(indXB[x2]))/nB)-gammaB[x1,y,z]/(max(1,length(indXB[x1]))/nB))
    @expression(modelB, regterm, sum(1/length(voisins_X[x1,:]) *reg_absB[x1,x2,y,z] for x1 in 1:nbX, x2 in findall(voisins_X[x1,:]), y in Y, z in Z))

    # by default, the OT cost and regularization term are weighted to lie in the same interval
    @objective(modelA, Min, sum(C[y,z] * gammaA[x,y,z]  for y in Y, z in Z, x in 1:nbX) + lambda_reg *  sum(1/length(voisins_X[x1,:]) *reg_absA[x1,x2,y,z] for x1 in 1:nbX, x2 in findall(voisins_X[x1,:]), y in Y, z in Z))

    @objective(modelB, Min, sum(C[y,z] * gammaB[x,y,z]  for y in Y, z in Z, x in 1:nbX) + lambda_reg *  sum(1/length(voisins_X[x1,:]) *reg_absB[x1,x2,y,z] for x1 in 1:nbX, x2 in findall(voisins_X[x1,:]), y in Y, z in Z))

   # Solve the problem
   optimize!(modelA)
   optimize!(modelB)

   # Extract the values of the solution
   gammaA_val = [value(gammaA[x,y,z]) for x in 1:nbX, y in Y, z in Z]
   gammaB_val = [value(gammaB[x,y,z]) for x in 1:nbX, y in Y, z in Z]

   # compute the resulting estimators for the distributions of Z conditional to X and Y in base A and of Y conditional to X and Z in base B
   estimatorZA = 1/length(Z) * ones(nbX,length(Y),length(Z))
   for x = 1:nbX
       for y in Y
           proba_c_mA = sum(gammaA_val[x,y,z] for z in Z)
           if proba_c_mA > 1.0e-6
                estimatorZA[x,y,:] = 1/proba_c_mA * gammaA_val[x,y,:]
            end
       end
   end
   estimatorYB = 1/length(Y) * ones(nbX,length(Y),length(Z))
   for x = 1:nbX
       for z in Z
           proba_c_mB = sum(gammaB_val[x,y,z] for y in Y)
           if proba_c_mB > 1.0e-6
                estimatorYB[x,:,z] = 1/proba_c_mB * gammaB_val[x,:,z]
            end
       end
   end

   # Display the solution
   # println("Solution of the joint probability transport")
   # println("Distance cost = ", sum(C[y,z] * (gammaA_val[x,y,z]+gammaB_val[x,y,z]) for y in Y, z in Z, x in 1:nbX))
   # println("Regularization cost = ", lambda_reg * value(regterm))


   Solution( time()-tstart,
             [sum(gammaA_val[x,y,z] for x in 1:nbX) for y in Y, z in Z],
             [sum(gammaB_val[x,y,z] for x in 1:nbX) for y in Y, z in Z],
              estimatorZA,
              estimatorYB)

end
