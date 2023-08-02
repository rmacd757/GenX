#### This file runs the initial paper runs for the inter-annual variability analysis

file_directory = @__DIR__

####### Run brownfields ########

include(joinpath(file_directory, "2zone_brownfields", "Run_fusion_1-5.jl"))

include(joinpath(file_directory, "2zone_brownfields", "Run_fusion_6-10.jl"))

include(joinpath(file_directory, "2zone_brownfields", "Run_fusion_11-15.jl"))

include(joinpath(file_directory, "2zone_brownfields", "Run_fusion_16-20.jl"))


include(joinpath(file_directory, "2zone_brownfields", "Run_nofusion_1-5.jl"))

include(joinpath(file_directory, "2zone_brownfields", "Run_nofusion_6-10.jl"))

include(joinpath(file_directory, "2zone_brownfields", "Run_nofusion_11-15.jl"))

include(joinpath(file_directory, "2zone_brownfields", "Run_nofusion_16-20.jl"))