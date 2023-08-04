#### This file runs the initial paper runs for the inter-annual variability analysis

file_directory = @__DIR__

####### Run 20 scenarios ########
include("/home/gridsan/nbhatt1/GenX/fusion_paper/paper_runs/2zone_20yr/fusion/Run.jl")

include("/home/gridsan/nbhatt1/GenX/fusion_paper/paper_runs/2zone_20yr/no_fusion/Run.jl")