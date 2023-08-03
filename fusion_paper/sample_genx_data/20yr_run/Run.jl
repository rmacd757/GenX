##

using GenX
using JuMP

case = dirname(@__FILE__)

function get_settings_path(case::AbstractString)
    return joinpath(case, "Settings")
end

function get_settings_path(case::AbstractString, filename::AbstractString)
    return joinpath(get_settings_path(case), filename)
end

function get_default_output_folder(case::AbstractString)
    return joinpath(case, "Results")
end

genx_settings = get_settings_path(case, "genx_settings.yml") #Settings YAML file path
mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters

inputs_path = case
settings_path = get_settings_path(case)

## Cluster time series inputs if necessary and if specified by the user
TDRpath = joinpath(inputs_path, mysetup["TimeDomainReductionFolder"])
if mysetup["TimeDomainReduction"] == 1
    if !time_domain_reduced_files_exist(TDRpath)
        println("Clustering Time Series Data (Grouped)...")
        cluster_inputs(inputs_path, settings_path, mysetup)
    else
        println("Time Series Data Already Clustered.")
    end
end

## Configure solver
println("Configuring Solver")
OPTIMIZER = configure_solver(mysetup["Solver"], settings_path)
# set_optimizer_attribute(OPTIMIZER, "BarHomogeneous", 1)

#### Running a case

## Load inputs
println("Loading Inputs")
myinputs = load_inputs(mysetup, inputs_path)

# emiss_lim_list = 100.0 .* [2.5, 5, 7.5, 10, 15, 20, 25]
# emiss_lim_list = 100.0 .* [2.5, 5, 10, 15, 20, 25, 1.0]
emiss_lim_list = 100.0 .* [5]

mysetup["CO2Cap"] = 1
scale_factor = mysetup["ParameterScale"] == 1 ? ModelScalingFactor : 1

for emiss_lim in emiss_lim_list

    # Hard-coded to put all emissions in New Hampshire, but CO2 Cap is set to be system-wide
    myinputs["dfMaxCO2"][2] = emiss_lim * 1e3 / scale_factor
    outputs_path = joinpath(case, "Results", "Primal_Capex_8500.0_EmissLevel_" * string(emiss_lim))

    if isfile(joinpath(outputs_path, "costs.csv"))
        println("Skipping Case for emiss limit = " * string(emiss_lim) * " because it already exists.")
        continue
    end

    mkpath(dirname(outputs_path))
    
    ## Generate model
    println("Generating the Optimization Model")
    EP = generate_model(mysetup, myinputs, OPTIMIZER)

    ########################
    #### Add any additional constraints
    HYDRO_RES = myinputs["HYDRO_RES"]
    dfGen = myinputs["dfGen"]

    ## Hydro storage <= 0.55 * Existing Capacity at start of May 1st 
    @constraint(EP, cHydroSpring[y in HYDRO_RES], EP[:vS_HYDRO][y, 2879] .<= 0.55 .* EP[:eTotalCap][y] .* dfGen[y,:Hydro_Energy_to_Power_Ratio]) 

    ## Hydro storage == 0.70 * Existing Capacity at the start of the year
    @constraint(EP, cHydroJan[y in HYDRO_RES], EP[:vS_HYDRO][y, 1] .== 0.70 .* EP[:eTotalCap][y] .* dfGen[y,:Hydro_Energy_to_Power_Ratio]) 

    ## Maine -> Quebec transmission limited to 2170MWe.
    # The line is defined as Quebec -> Maine in Network.csv, so these flows will be negative
    # Make sure to correc the line index if the order is changed in Network.csv
    @constraint(EP, cMaine2Quebec[t=1:myinputs["T"]], EP[:vFLOW][2, t] >= -170.0)

    ## Solar <= 22GWe
    # solar_rid = findall(x -> startswith(x, "solar"), dfGen[!,:Resource])
    # @constraint(EP, cSolarCap, sum(EP[:eTotalCap][y] for y in solar_rid) <= 22e3)

    # ## Onshore wind <= 10GWe
    # onshore_rid = findall(x -> startswith(x, "onshore"), dfGen[!,:Resource])
    # @constraint(EP, cOnshoreCap, sum(EP[:eTotalCap][y] for y in onshore_rid) <= 10e3)

    # ## Offshore wind <= 280GWe
    # offshore_rid = findall(x -> startswith(x, "offshore"), dfGen[!,:Resource])
    # @constraint(EP, cOffshoreCap, sum(EP[:eTotalCap][y] for y in offshore_rid) <= 280e3)

    ## Offshore_fixed <= 37.5 GWe 
    offshore_fixed = findall(x -> startswith(x, "fixed_offshore"), dfGen[!,:Resource])
    @constraint(EP, cOffshoreFixCap, sum(EP[:eTotalCap][y] for y in offshore_fixed) <= 37.5e3)

    ## Offshore_floating <= 275 GWe 
    offshore_float = findall(x -> startswith(x, "float_offshore"), dfGen[!,:Resource])
    @constraint(EP, cOffshoreFloatCap, sum(EP[:eTotalCap][y] for y in offshore_float) <= 275e3)

    ########################

    ## Solve model
    println("Solving Model")
    EP, solve_time = solve_model(EP, mysetup)
    myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

    ## Run MGA if the MGA flag is set to 1 else only save the least cost solution
    println("Writing Output")
    # outputs_path = get_default_output_folder(outputs_path)

    ## Write outputs
    write_outputs(EP, outputs_path, mysetup, myinputs)

end






















# ### Cluster time series inputs if necessary and if specified by the user
# TDRpath = joinpath(case, mysetup["TimeDomainReductionFolder"])

# function time_domain_reduced_files_exist(tdrpath)
#     tdr_load = isfile(joinpath(tdrpath,"Load_data.csv"))
#     tdr_genvar = isfile(joinpath(tdrpath,"Generators_variability.csv"))
#     tdr_fuels = isfile(joinpath(tdrpath,"Fuels_data.csv"))
#     return (tdr_load && tdr_genvar && tdr_fuels)
# end

# if mysetup["TimeDomainReduction"] == 1
#     if !time_domain_reduced_files_exist(TDRpath)
#         println("Clustering Time Series Data (Grouped)...")
#         cluster_inputs(inputs_path, settings_path, mysetup)
#     else
#         println("Time Series Data Already Clustered.")
#     end
# end

# ### Configure solver
# println("Configuring Solver")
# OPTIMIZER = configure_solver(mysetup["Solver"], settings_path)


# #### Running a case

# ### Load inputs
# println("Loading Inputs")
# myinputs = load_inputs(mysetup, case)

# println("Generating the Optimization Model")
# EP = generate_model(mysetup, myinputs, OPTIMIZER)

# ########################
# #### Add any additional constraints
# HYDRO_RES = myinputs["HYDRO_RES"]
# dfGen = myinputs["dfGen"]

# ## Hydro storage <= 0.55 * Existing Capacity at start of May 1st 
# @constraint(EP, cHydroSpring[y in HYDRO_RES], EP[:vS_HYDRO][y, 2879] .<= 0.55 .* EP[:eTotalCap][y] .* dfGen[y,:Hydro_Energy_to_Power_Ratio]) 

# ## Hydro storage == 0.70 * Existing Capacity at the start of the year
# @constraint(EP, cHydroJan[y in HYDRO_RES], EP[:vS_HYDRO][y, 1]       .== 0.70 .* EP[:eTotalCap][y] .* dfGen[y,:Hydro_Energy_to_Power_Ratio]) 

# ## Maine -> Quebec transmission limited to 2170MWe.
# # The line is defined as Quebec -> Maine in Network.csv, so these flows will be negative
# # Make sure to correc the line index if the order is changed in Network.csv
# @constraint(EP, cMaine2Quebec[t=1:myinputs["T"]], EP[:vFLOW][2, t] >= -170.0)

# ## Solar <= 22GWe
# solar_rid = findall(x -> startswith(x, "solar"), dfGen[!,:Resource])
# @constraint(EP, cSolarCap, sum(EP[:eTotalCap][y] for y in solar_rid) <= 22e3)

# ## Onshore wind <= 10GWe
# onshore_rid = findall(x -> startswith(x, "onshore"), dfGen[!,:Resource])
# @constraint(EP, cOnshoreCap, sum(EP[:eTotalCap][y] for y in onshore_rid) <= 10e3)

# ## Offshore wind <= 280GWe
# offshore_rid = findall(x -> startswith(x, "offshore"), dfGen[!,:Resource])
# @constraint(EP, cOffshoreCap, sum(EP[:eTotalCap][y] for y in offshore_rid) <= 280e3)

# ########################


# println("Solving Model")
# EP, solve_time = solve_model(EP, mysetup)
# myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

# # Run MGA if the MGA flag is set to 1 else only save the least cost solution
# println("Writing Output")
# outputs_path = get_default_output_folder(case)
# elapsed_time = @elapsed write_outputs(EP, outputs_path, mysetup, myinputs)
# println("Time elapsed for writing is")
# println(elapsed_time)

