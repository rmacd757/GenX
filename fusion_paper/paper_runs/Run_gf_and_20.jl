#### This file runs the initial paper runs for the inter-annual variability analysis

file_directory = @__DIR__

####### Run greenfields ########
include(joinpath(file_directory, "2zone_greenfields", "Run_nofus.jl"))

include(joinpath(file_directory, "2zone_greenfields", "Run_fus.jl"))

####### Run 20 scenarios ########
include(joinpath(file_directory, "2zone_20yr/89QC_Flex_Fus", "Run.jl"))

inlcude(joinpath(file_directory, "2zone_20yr/89QC_Flex_NoFus", "Run.jl"))