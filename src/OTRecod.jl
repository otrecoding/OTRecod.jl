module OTRecod

using Statistics
using JuMP, Cbc, Clp
using Printf
using DataFrames
using Distances
using DocStringExtensions

export run_directory
export ot_group, ot_joint

include("instance.jl")
include("solution.jl")
include("simulate.jl")
include("disp_inst_inf.jl")
include("empirical_distribution.jl")
include("bound_prediction_error.jl")
include("aggregate_per_covar_mixed.jl")
include("average_distance_closest.jl")
include("distrib_error.jl")
include("pred_error.jl")
include("ot_group.jl")
include("ot_joint.jl")
include("compute_average_error_bound.jl")
include("empirical_estimator.jl")
include("run_directory.jl")

end
