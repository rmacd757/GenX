using GenX
using JuMP
using OrderedCollections
using DataFrames
using CSV

input_name = "final_thermstor_singleyear"
case_name = "final_thermstor_singleyear"

function gethomedir(case_path::String)
    path_split = splitpath(case_path)
    home_dir = ""
    for s in path_split
        if s == "fusion_paper"
            home_dir = joinpath(home_dir, s)
            break
        end
        home_dir = joinpath(home_dir, s)
    end

    return home_dir
end

function set_all_thermstor_prop!(dfFusion::DataFrame, fusion_rid::Vector{Int}, thermstor_type::String, location_adjustment::Float64)
    turb_cost = 1700.0
    for y in fusion_rid
        turb_cost = set_thermstor_prop!(dfFusion, y, thermstor_type, location_adjustment)
    end
    return turb_cost
end

function set_thermstor_prop!(dfFusion::DataFrame, y::Int, thermstor_type::String, location_adjustment::Float64)
    if thermstor_type == "final_thermstor"
        dfFusion[y,:Stor_Cost_per_MWht] = 891.9610585 
        dfFusion[y,:Dis_Cost_per_MWht] = 7433.00882
        dfFusion[y,:Therm_Stor_Leakage_Daily_Frac] = 0.02
        turb_cost = 1700.0
        dfFusion[y,:Turb_CAPEX] = turb_cost .* 1000.0 .* location_adjustment
    elseif thermstor_type == "final_thermstor_noLeak"
        dfFusion[y,:Stor_Cost_per_MWht] = 891.9610585 
        dfFusion[y,:Dis_Cost_per_MWht] = 7433.00882
        dfFusion[y,:Therm_Stor_Leakage_Daily_Frac] = 0.0
        turb_cost = 1700.0
        dfFusion[y,:Turb_CAPEX] = turb_cost .* 1000.0 .* location_adjustment
    elseif thermstor_type == "final_thermstor_noLeak_cheapTurb"
        dfFusion[y,:Stor_Cost_per_MWht] = 891.9610585 
        dfFusion[y,:Dis_Cost_per_MWht] = 7433.00882
        dfFusion[y,:Therm_Stor_Leakage_Daily_Frac] = 0.0
        turb_cost = 1000.0
        dfFusion[y,:Turb_CAPEX] = turb_cost .* 1000.0 .* location_adjustment
    elseif thermstor_type == "final_thermstor_noLeak_cheapTurb750"
        dfFusion[y,:Stor_Cost_per_MWht] = 891.9610585 
        dfFusion[y,:Dis_Cost_per_MWht] = 7433.00882
        dfFusion[y,:Therm_Stor_Leakage_Daily_Frac] = 0.0
        turb_cost = 750.0
        dfFusion[y,:Turb_CAPEX] = turb_cost .* 1000.0 .* location_adjustment
    elseif thermstor_type == "final_thermstor_10"
        dfFusion[y,:Stor_Cost_per_MWht] = 89.19610585 
        dfFusion[y,:Dis_Cost_per_MWht] = 743.300882
        dfFusion[y,:Therm_Stor_Leakage_Daily_Frac] = 0.02
        turb_cost = 1700.0
        dfFusion[y,:Turb_CAPEX] = turb_cost .* 1000.0 .* location_adjustment
    elseif thermstor_type == "final_thermstor_10_noLeak"
        dfFusion[y,:Stor_Cost_per_MWht] = 89.19610585 
        dfFusion[y,:Dis_Cost_per_MWht] = 743.300882
        dfFusion[y,:Therm_Stor_Leakage_Daily_Frac] = 0.0
        turb_cost = 1700.0
        dfFusion[y,:Turb_CAPEX] = turb_cost .* 1000.0 .* location_adjustment
    elseif thermstor_type == "final_thermstor_10_noLeak_cheapTurb"
        dfFusion[y,:Stor_Cost_per_MWht] = 89.19610585 
        dfFusion[y,:Dis_Cost_per_MWht] = 743.300882
        dfFusion[y,:Therm_Stor_Leakage_Daily_Frac] = 0.0
        turb_cost = 1000.0
        dfFusion[y,:Turb_CAPEX] = turb_cost .* 1000.0 .* location_adjustment
    elseif thermstor_type == "final_thermstor_10_noLeak_cheapTurb750"
        dfFusion[y,:Stor_Cost_per_MWht] = 89.19610585 
        dfFusion[y,:Dis_Cost_per_MWht] = 743.300882
        dfFusion[y,:Therm_Stor_Leakage_Daily_Frac] = 0.0
        turb_cost = 750.0
        dfFusion[y,:Turb_CAPEX] = turb_cost .* 1000.0 .* location_adjustment
    else
        println("Thermstor type not recognized")
        turb_cost = 1700.0
    end
    return turb_cost
end

# Find the home directory, to let us load the run_helpers.jl file
case_path = @__DIR__
home_dir = gethomedir(case_path)
println(home_dir) 

## Load helper functions
include(joinpath(home_dir,"run_helpers.jl"))

# Total load = 4,827,887,023 MWh across all 20 scenarios
# The scenarios vary in load by ~1%, so we'll treat the equally
# The emission intensity limits are: [4.0, 12.0, 50.0] gCO2 / kWh
# GenX requires the limit to be in millions tonnes (metric), so we'll convert by:
# g / tonne = 1e6
# kWh / MWh = 1e3
# total = 4,827,887,023[MWh] * limit[g/kWh] * 1e3[kWh/MWh] / 1e6[tonne/g] / 1e6[MMT / tonne]
# total = 4,827,887,023 * limit / 1e9
# emiss_lim_list = [4.0, 12.0, 50.0]
emiss_lim_list = [0.0, 2.5]

fusion_cost_list = [8500.0, 3000.0, 6000.0, 12000.0]

thermstor_type_list = [
    "final_thermstor",
    "final_thermstor_noLeak",
    "final_thermstor_noLeak_cheapTurb",
    "final_thermstor_noLeak_cheapTurb750",
    "final_thermstor_10",
    "final_thermstor_10_noLeak",
    "final_thermstor_10_noLeak_cheapTurb",
    "final_thermstor_10_noLeak_cheapTurb750"
    
]

scenario_list = 0:1:19

task_id = parse(Int,ARGS[1])
num_tasks = parse(Int,ARGS[2])
num_threads = parse(Int,ARGS[3])

all_cases = vcat(collect(Iterators.product(emiss_lim_list, fusion_cost_list, thermstor_type_list, scenario_list))...)

reduced_cases = []

# Go through the cases and add any where !isfile(joinpath(outputs_path, "costs.csv"))
# for idx in task_id+1:num_tasks:length(all_cases)
#     emiss_lim = all_cases[idx][1]
#     fusion_cost = all_cases[idx][2]
#     thermstor_type = all_cases[idx][3]
#     scenario_number = all_cases[idx][4]

#     mkpath(joinpath(results_path, thermstor_type))

#     outputs_path = joinpath(results_path, thermstor_type, "Cost_$(fusion_cost)_EmissLevel_$(emiss_lim)_gCO2perkWh_Scenario_$(scenario_number)")
#     if !isfile(joinpath(outputs_path, "costs.csv"))
#         println("Including Case for emiss limit = $emiss_lim, fusion cap = $fusion_cost, thermstor type = $thermstor_type, scenario number = $scenario_number")
#         push!(reduced_cases, (emiss_lim, fusion_cost, thermstor_type, scenario_number))
#         rm(outputs_path, force=true, recursive=true)
#     end
# end

for idx in task_id+1:num_tasks:length(all_cases)
    GC.gc()

    emiss_lim = all_cases[idx][1]
    fusion_cost = all_cases[idx][2]
    thermstor_type = all_cases[idx][3]
    scenario_number = all_cases[idx][4]

    ## Define input and output paths
    inputs_path = joinpath(@__DIR__, "Scenarios", "scenario_$(scenario_number)")
    results_path = joinpath(@__DIR__, "Results")

    # mkpath(results_path)
    mkpath(joinpath(results_path, thermstor_type))

    ## Load settings
    genx_settings = get_settings_path(inputs_path, "genx_settings.yml") #Settings YAML file path
    mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters
    settings_path = get_settings_path(inputs_path)

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

    mysetup["CO2Cap"] = 2
    scale_factor = mysetup["ParameterScale"] == 1 ? ModelScalingFactor : 1

    ## Configure solver
    println("Configuring Solver")
    OPTIMIZER = configure_solver(mysetup["Solver"], settings_path)
    set_optimizer_attribute(OPTIMIZER, "Threads", num_threads)

    # Turn this setting on if you run into numerical stability issues
    # set_optimizer_attribute(OPTIMIZER, "BarHomogeneous", 1)

    #### Running a case

    ## Load inputs
    println("Loading Inputs")
    myinputs = load_inputs(mysetup, inputs_path)

    println("Emiss Limit: $emiss_lim, Fusion Cost: $fusion_cost, Thermstor Type: $thermstor_type, Scenario Number: $scenario_number")

    myinputs["dfMaxCO2Rate"][2] = emiss_lim / scale_factor ./ 1e3
    outputs_path = joinpath(results_path, thermstor_type, "Cost_$(fusion_cost)_EmissLevel_$(emiss_lim)_gCO2perkWh_Scenario_$(scenario_number)")
    mkpath(outputs_path)

    discount_factor = 0.06
    lifetime = 40.0
    annuity = discount_factor / (1.0 - (1.0 + discount_factor)^(-lifetime))
    # turb_cost = 1700.0
    vessel_cost = 150.0
    location_adjustment = 1.12
    fixed_cost_ratio = 0.15
    num_years = 1.0

    dfFusion = myinputs["dfFusion"]
    fusion_rid = findall(x -> startswith(x, "fusion"), dfFusion[!,:Resource])
    turb_cost = set_all_thermstor_prop!(dfFusion, fusion_rid, thermstor_type, location_adjustment)

    fusion_annual_cost = (fusion_cost .- turb_cost .- vessel_cost) .* annuity .* num_years .* 1000 .* location_adjustment
    fusion_fixed_cost = fusion_cost .* annuity .* num_years .* 1000 .* fixed_cost_ratio

    # Find all the fusion resources in the model
    # and set their investment and fixed O&M costs to zero
    dfGen = myinputs["dfGen"]
    fusion_rid = findall(x -> startswith(x, "fusion"), dfGen[!,:Resource])
    for y in fusion_rid
        println("Running FPP at R_ID $y with investment costs: $(fusion_annual_cost) and fixed costs: $(fusion_fixed_cost)")
        dfGen[y,:Inv_Cost_per_MWyr] = fusion_annual_cost 
        dfGen[y,:Fixed_OM_Cost_per_MWyr] = fusion_fixed_cost
    end

    # This check will cause the case to be skipped if the results already exist
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

    # Empty arrays for indexing
    jan1_idxs = Int[]
    may1_idxs = Int[]

    # 20 year indexing
    for year_num in 1:num_years
        # Calculate the index for the beginning of the years
        start_year = (year_num-1) * 8760 + 1
        push!(jan1_idxs, start_year)

        # Calculate the index for the middle of the years
        mid_year = (year_num-1) * 8760 + 2879
        push!(may1_idxs, mid_year)
    end

    ## Hydro storage == 0.70 * Existing Capacity at the start of the year
    @constraint(EP, cHydroJan[y in HYDRO_RES, jan1_idx in jan1_idxs], EP[:vS_HYDRO][y, jan1_idx]  .== 0.70 .* EP[:eTotalCap][y] .* dfGen[y,:Hydro_Energy_to_Power_Ratio])
    
    ## Hydro storage <= 0.55 * Existing Capacity at start of May 1st 
    @constraint(EP, cHydroSpring[y in HYDRO_RES, may1_idx in may1_idxs], EP[:vS_HYDRO][y, may1_idx] .<= 0.55 .* EP[:eTotalCap][y] .* dfGen[y,:Hydro_Energy_to_Power_Ratio])
        
    ## Maine -> Quebec transmission limited to 2170MWe.
    # The line is defined as Quebec -> Maine in Network.csv, so these flows will be negative
    # Make sure to correc the line index if the order is changed in Network.csv
    @constraint(EP, cMaine2Quebec[t=1:myinputs["T"]], EP[:vFLOW][2, t] >= -170.0)

    ########################

    ## Solve model
    println("Solving Model")
    EP, solve_time = solve_model(EP, mysetup)
    myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

    ## Run MGA if the MGA flag is set to 1 else only save the least cost solution
    println("Writing Output")
    # outputs_path = get_default_output_folder(outputs_path)

    ## Write outputs
    # write_outputs(EP, outputs_path, mysetup, myinputs)
    GenX.write_costs(outputs_path, myinputs, mysetup, EP)
    GenX.write_capacity(outputs_path, myinputs, mysetup, EP)
    GenX.write_capacityfactor(outputs_path, myinputs, mysetup, EP)
    # GenX.write_nse(outputs_path, myinputs, mysetup, EP)
    GenX.write_fusion_summary(outputs_path, myinputs, mysetup, EP)
    # result_summ = DataFrame(Cost=objective_value(EP), Dual=0.0)
    # CSV.write(joinpath(outputs_path, "fpp_results.csv"), result_summ)

end

