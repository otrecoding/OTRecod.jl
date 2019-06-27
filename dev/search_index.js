var documenterSearchIndex = {"docs":
[{"location":"OT_group/#Group-1","page":"Group","title":"Group","text":"","category":"section"},{"location":"OT_group/#","page":"Group","title":"Group","text":"Modules = [OTRecod]\nPages   = [\"OT_group.jl\"]","category":"page"},{"location":"OT_group/#OTRecod.OT_group","page":"Group","title":"OTRecod.OT_group","text":"OT_group(inst, percent_closest=0.2, maxrelax=0.0, norme=0, \n         indiv_method=sequential, full_disp=false, solver_disp=false)\n\nModel of group transport percent_closest: percent of closest neighbors taken in the computation of the   costs\n\nmaxrelax: maximum percentage of deviation from expected probability masses\nindiv_method: specifies the method used to get individual transport from group joint probabilities\nfull_disp: if true, write the transported value of each individual;     otherwise, juste write the number of missed transports\nsolver_disp: if false, do not display the outputs of the solver\n\n\n\n\n\n","category":"function"},{"location":"OT_group/#OTRecod.individual_from_group_closest","page":"Group","title":"OTRecod.individual_from_group_closest","text":"individual_from_group_closest(inst, jointprobaA, \n    jointprobaB, percent_closest=1.0)\n\nSequentially assign the modality of the individuals to that of the closest neighbor in the other base until the joint probability values are met\n\n\n\n\n\n","category":"function"},{"location":"OT_group/#OTRecod.individual_from_group_optimal","page":"Group","title":"OTRecod.individual_from_group_optimal","text":"individual_from_group_optimal(inst, jointprobaA,\n     jointprobaB; percent_closest=1.0)\n\nSolve an optimization problem to get the individual transport that minimizes total distance while satisfying the joint probability computed by the model by group\n\n\n\n\n\n","category":"function"},{"location":"OT_joint/#Joint-1","page":"Joint","title":"Joint","text":"","category":"section"},{"location":"OT_joint/#","page":"Joint","title":"Joint","text":"Modules = [OTRecod]\nPages   = [\"OT_joint.jl\"]","category":"page"},{"location":"OT_joint/#OTRecod.OT_joint","page":"Joint","title":"OTRecod.OT_joint","text":"Model where we directly compute the distribution of the outcomes for each individual or for sets of indviduals that similar values of covariates aggregatetol: quantify how much individuals' covariates must be close for aggregation regnorm: norm1, norm2 or entropy depending on the type of regularization percentclosest: percent of closest neighbors taken into consideration in regularization lambdareg: coefficient measuing the importance of the regularization term fulldisp: if true, write the transported value of each individual; otherwise, juste write the number of missed transports solverdisp: if false, do not display the outputs of the solver\n\n\n\n\n\n","category":"function"},{"location":"PrintLog/#Logging-1","page":"Logging","title":"Logging","text":"","category":"section"},{"location":"PrintLog/#","page":"Logging","title":"Logging","text":"Modules = [OTRecod]\nPages   = [\"PrintLog.jl\"]","category":"page"},{"location":"#OTRecod.jl-1","page":"OTRecod.jl","title":"OTRecod.jl","text":"","category":"section"},{"location":"#","page":"OTRecod.jl","title":"OTRecod.jl","text":"Documentation for OTRecod.jl","category":"page"},{"location":"plot_functions/#Plots-1","page":"Plots","title":"Plots","text":"","category":"section"},{"location":"plot_functions/#","page":"Plots","title":"Plots","text":"Modules = [OTRecod]\nPages   = [\"plot_functions.jl\"]","category":"page"},{"location":"utils/#Helper-functions-1","page":"Helper functions","title":"Helper functions","text":"","category":"section"},{"location":"utils/#","page":"Helper functions","title":"Helper functions","text":"Modules = [OTRecod]\nPages   = [\"utils.jl\"]","category":"page"},{"location":"utils/#OTRecod.Instance","page":"Helper functions","title":"OTRecod.Instance","text":"Definition and initialization of an Instance structure\n\n\n\n\n\n","category":"type"},{"location":"utils/#OTRecod.average_distance_to_closest-Tuple{Instance,Float64}","page":"Helper functions","title":"OTRecod.average_distance_to_closest","text":"Compute the cost between pairs of outcomes as the average distance between covariations of individuals with these outcomes, but considering only the percent closest neighbors\n\n\n\n\n\n","category":"method"},{"location":"utils/#OTRecod.avg_distance_closest-Tuple{Instance,OTRecod.DataBase,OTRecod.DataBase,OTRecod.DataBase,Int64,Int64,Float64}","page":"Helper functions","title":"OTRecod.avg_distance_closest","text":"Compute the average distance between individuals of base1 with modality m1 for outcome and individuals of base2 with modality m2 for outcome Consider only the percent_closest individuals in the computation of the distance\n\n\n\n\n\n","category":"method"},{"location":"utils/#OTRecod.bound_prediction_error","page":"Helper functions","title":"OTRecod.bound_prediction_error","text":"Compute a bound on the average prediction error in each base The bound is computed as the expected prediction error assuming that the distribution of Z in base A (and that of Y in base B) is known, and the prediction done with the value that maximizes the probability\n\n\n\n\n\n","category":"function"},{"location":"utils/#OTRecod.compute_pred_error","page":"Helper functions","title":"OTRecod.compute_pred_error","text":"Compute prediction errors in a solution\n\n\n\n\n\n","category":"function"},{"location":"utils/#OTRecod.disp_inst_info-Tuple{Instance}","page":"Helper functions","title":"OTRecod.disp_inst_info","text":"Display information about the distance between the modalities\n\n\n\n\n\n","category":"method"},{"location":"utils/#OTRecod.empirical_distribution","page":"Helper functions","title":"OTRecod.empirical_distribution","text":"Return the empirical cardinality of the joint occurrences of (C=x,Y=mA,Z=mB) in both bases\n\n\n\n\n\n","category":"function"}]
}
