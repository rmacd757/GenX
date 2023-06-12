home_dir = dirname(@__FILE__)
include(joinpath(home_dir, "primal_2zone_20year_flexNG_fppThermStor_highcost", "Run.jl"))
include(joinpath(home_dir, "primal_2zone_20year_flexNG_baseloadfpp", "Run.jl"))
