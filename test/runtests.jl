using OTRecod, Test
using MultivariateStats
using Distributions, PDMats
using Printf

include("test_utils.jl")
include("test_params_group.jl")
include("test_params_joint.jl")
#include("test_ncds.jl")

@testset "run directory " begin

    for method in [:group, :joint]

        outfile = "result_$(method).out"

        tcpu = @time run_directory("data", method, outfile)

        open(outfile,"r") do f
            for line in eachline(f)
                print(line)
            end
        end

        @show tcpu
        @test true

    end

end
