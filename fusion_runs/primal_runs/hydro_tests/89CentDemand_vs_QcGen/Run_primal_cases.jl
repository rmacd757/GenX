home_dir = dirname(@__FILE__)
include(joinpath(home_dir, "hydrotest_89QcDemand_primal_6zoneandQC_noEmissLim", "Run.jl"))
include(joinpath(home_dir, "hydrotest_89QcDemand_primal_6zoneandQC_baseline", "Run.jl"))
include(joinpath(home_dir, "hydrotest_89QcDemand_primal_6zoneandQC", "Run.jl"))
