var documenterSearchIndex = {"docs":
[{"location":"#OTRecod.jl-1","page":"Documentation","title":"OTRecod.jl","text":"","category":"section"},{"location":"#","page":"Documentation","title":"Documentation","text":"Documentation for OTRecod.jl","category":"page"},{"location":"#Installation-1","page":"Documentation","title":"Installation","text":"","category":"section"},{"location":"#","page":"Documentation","title":"Documentation","text":"In a Julia session switch to pkg> mode to add NPSMC:","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"julia>] # switch to pkg> mode\npkg> add https://github.com/otrecoding/OTRecod.jl","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"Alternatively, you can achieve the above using the Pkg API:","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"julia> using Pkg\njulia> pkg\"add https://github.com/otrecoding/OTRecod.jl\"","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"When finished, make sure that you're back to the Julian prompt (julia>) and bring OTRecod into scope:","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"julia> using OTRecod","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"You can test the package with","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"julia>] # switch to pkg> mode\npkg> test OTRecod","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"To run an example from a dataset","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"run_directory","category":"page"},{"location":"#OTRecod.run_directory","page":"Documentation","title":"OTRecod.run_directory","text":"run_directory(path, method; outname=\"result.out\", \n                            maxrelax=0.0, \n                            lambda_reg=0.0, \n                            nbfiles=0, \n                            norme=0, \n                            percent_closest=0.2)\n\nRun one given method on a given number of data files of a given directory  The data files must be the only files with extension \".txt\" in the directory\n\npath   : name of the directory\nmethod : :group or :joint\nnbfiles: number of files considered, 0 if all the data files are tested\nnorme  : 1 or 2, norm used for distances in the space of covariates\n\n(see runallmethods for the description of other parameters)\n\n\n\n\n\n","category":"function"},{"location":"#","page":"Documentation","title":"Documentation","text":"julia> using OTRecod","category":"page"},{"location":"OT_group/#Group-1","page":"Group","title":"Group","text":"","category":"section"},{"location":"OT_group/#","page":"Group","title":"Group","text":"Modules = [OTRecod]\nPages   = [\"OT_group.jl\"]","category":"page"},{"location":"OT_group/#OTRecod.OT_group","page":"Group","title":"OTRecod.OT_group","text":"OT_group(inst, percent_closest=0.2, maxrelax=0.0, norme=0, \n         indiv_method=sequential, full_disp=false, solver_disp=false)\n\nModel of group transport percent_closest: percent of closest neighbors taken in the computation of the   costs\n\nmaxrelax: maximum percentage of deviation from expected probability masses\nindiv_method: specifies the method used to get individual transport from group joint probabilities\nfull_disp: if true, write the transported value of each individual;     otherwise, juste write the number of missed transports\nsolver_disp: if false, do not display the outputs of the solver\n\n\n\n\n\n","category":"function"},{"location":"OT_group/#OTRecod.individual_from_group_closest","page":"Group","title":"OTRecod.individual_from_group_closest","text":"individual_from_group_closest(inst, jointprobaA, jointprobaB, percent_closest=1.0)\n\nSequentially assign the modality of the individuals to that of the closest neighbor in the other base until the joint probability values are met\n\n\n\n\n\n","category":"function"},{"location":"OT_group/#OTRecod.individual_from_group_optimal","page":"Group","title":"OTRecod.individual_from_group_optimal","text":"individual_from_group_optimal(inst, jointprobaA,\n     jointprobaB; percent_closest=1.0)\n\nSolve an optimization problem to get the individual transport that minimizes total distance while satisfying the joint probability computed by the model by group\n\n\n\n\n\n","category":"function"},{"location":"OT_joint/#Joint-1","page":"Joint","title":"Joint","text":"","category":"section"},{"location":"OT_joint/#","page":"Joint","title":"Joint","text":"Modules = [OTRecod]\nPages   = [\"OT_joint.jl\"]","category":"page"},{"location":"OT_joint/#OTRecod.OT_joint","page":"Joint","title":"OTRecod.OT_joint","text":"Model where we directly compute the distribution of the outcomes for each individual or for sets of indviduals that similar values of covariates\n\naggregate_tol: quantify how much individuals' covariates must be close for aggregation\nreg_norm: norm1, norm2 or entropy depending on the type of regularization\npercent_closest: percent of closest neighbors taken into consideration in regularization\nlambda_reg: coefficient measuing the importance of the regularization term\nfull_disp: if true, write the transported value of each individual; otherwise, juste write the number of missed transports\nsolver_disp: if false, do not display the outputs of the solver\n\n\n\n\n\n","category":"function"},{"location":"PrintLog/#Logging-1","page":"Logging","title":"Logging","text":"","category":"section"},{"location":"PrintLog/#","page":"Logging","title":"Logging","text":"Modules = [OTRecod]\nPages   = [\"PrintLog.jl\"]","category":"page"},{"location":"plot_functions/#Plots-1","page":"Plotting","title":"Plots","text":"","category":"section"},{"location":"plot_functions/#","page":"Plotting","title":"Plotting","text":"Modules = [OTRecod]\nPages   = [\"plot_functions.jl\"]","category":"page"},{"location":"plot_functions/#OTRecod.plot_benchmark","page":"Plotting","title":"OTRecod.plot_benchmark","text":"plot_benchmark(datapath,outputpath,ymin=0.0, ymax=0.8)\n\nBoxplot of each simulation with all methods\n\n\n\n\n\n","category":"function"},{"location":"utils/#Helper-functions-1","page":"Utilities","title":"Helper functions","text":"","category":"section"},{"location":"utils/#","page":"Utilities","title":"Utilities","text":"Modules = [OTRecod]\nPages   = [\"utils.jl\"]","category":"page"},{"location":"utils/#OTRecod.average_distance_to_closest-Tuple{Instance,Float64}","page":"Utilities","title":"OTRecod.average_distance_to_closest","text":"average_distance_to_closest(instance, percent_closest)\n\nCompute the cost between pairs of outcomes as the average distance between covariations of individuals with these outcomes, but considering only the percent closest neighbors\n\n\n\n\n\n","category":"method"},{"location":"utils/#OTRecod.avg_distance_closest-Tuple{Instance,OTRecod.DataBase,OTRecod.DataBase,OTRecod.DataBase,Int64,Int64,Float64}","page":"Utilities","title":"OTRecod.avg_distance_closest","text":"avg_distance_closest(instance, database1, database2, outcome, m1, m2, \n                     percent_closest)\n\nCompute the average distance between individuals of base1 with modality m1 for outcome and individuals of base2 with modality m2 for outcome Consider only the percent_closest individuals in the computation of the distance\n\n\n\n\n\n","category":"method"},{"location":"utils/#OTRecod.bound_prediction_error","page":"Utilities","title":"OTRecod.bound_prediction_error","text":"bound_prediction_error(instance, norme, aggregate_tol)\n\nCompute a bound on the average prediction error in each base. The bound is computed as the expected prediction error assuming that the distribution of Z in base A (and that of Y in base B) is known, and the prediction done with the value that maximizes the probability\n\n\n\n\n\n","category":"function"},{"location":"utils/#OTRecod.compute_distrib_error-Tuple{Instance,OTRecod.Solution,Any,Any}","page":"Utilities","title":"OTRecod.compute_distrib_error","text":"compute_distrib_error(insttance, solution, empiricalZA, empiricalYB)\n\nCompute errors in the conditional distributions of a solution\n\n\n\n\n\n","category":"method"},{"location":"utils/#OTRecod.compute_pred_error","page":"Utilities","title":"OTRecod.compute_pred_error","text":"Compute prediction errors in a solution\n\n\n\n\n\n","category":"function"},{"location":"utils/#OTRecod.disp_inst_info-Tuple{Instance}","page":"Utilities","title":"OTRecod.disp_inst_info","text":"disp_inst_info(instance)\n\nDisplay information about the distance between the modalities\n\n\n\n\n\n","category":"method"},{"location":"utils/#OTRecod.empirical_distribution","page":"Utilities","title":"OTRecod.empirical_distribution","text":"empirical_distribution(instance, norme, aggregate_tol)\n\nReturn the empirical cardinality of the joint occurrences of (C=x,Y=mA,Z=mB) in both bases\n\n\n\n\n\n","category":"function"}]
}
