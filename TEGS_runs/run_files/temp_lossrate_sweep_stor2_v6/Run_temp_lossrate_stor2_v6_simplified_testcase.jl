using GenX
using CSV
using DataFrames
using JuMP

# include(joinpath(pwd(),"PiecewiseLinearOpt","PiecewiseLinearOpt.jl"))

function tegs_case(dsratio, csratio, max_temp, min_temp, tegs_rid, outputs_path_case, mysetup, myinputs, OPTIMIZER)
    dfGen = myinputs["dfGen"]

    START_SUBPERIODS = myinputs["START_SUBPERIODS"]
    INTERIOR_SUBPERIODS = myinputs["INTERIOR_SUBPERIODS"]
    hours_per_subperiod = myinputs["hours_per_subperiod"]

    # Calculate and save baseline emissions
    mkpath(outputs_path_case)
    outputs_path = joinpath(outputs_path_case, "ds-$(dsratio)_cs-$(dsratio)")

    # if isdir(outputs_path)
    #     println("Result already exists")
    #     temp = CSV.File(joinpath(outputs_path, "costs.csv"))::CSV.File
    #     return temp[1][2]
    # end
    
    println("Generating the Optimization Model")
    EP = generate_model(mysetup, myinputs, OPTIMIZER)

    @constraint(EP, min_tegs, EP[:eTotalCapEnergy][end] >= 10000.0)
    # @constraint(EP, start_vS, EP[:vS][tegs_rid,1] == -100.0)

    ## Fix discharge capacity - storage capacity ratio
    # {Discharge capacity} * {dsratio} == {energy capacity}
    @constraint(EP, TEGSdsratio, EP[:eTotalCap][end] * dsratio == EP[:eTotalCapEnergy][end])

    ## Fix charging capacity - storage capacity ratio
    # {Charge capacity} * {csratio} == {energy capacity}
    @constraint(EP, TEGScsratio, EP[:eTotalCapCharge][end] * csratio == EP[:eTotalCapEnergy][end])

    ## Allow TEGS to go to a negative charge state
    for t in 1:myinputs["T"]
        delete_lower_bound(EP[:vS][end,t])
    end
    # Remove existing energy balance constraints
    for t in START_SUBPERIODS
        delete(EP, EP[:cSoCBalStart][t, end])
        delete(EP, EP[:cStart_vP_leq_vS][end,t])
        # unregister(EP, EP[:cSoCBalStart][t, end])
    end
    for t in INTERIOR_SUBPERIODS
        delete(EP, EP[:cSoCBalInterior][t, end])
        delete(EP, EP[:cInterior_vP_leq_vS][end,t])
        # unregister(EP, EP[:cSoCBalInterior][t, end])
    end

    ## Make leakage continue when state of charge = 0% but temperature > ambient temperature
    AMBIENT_TEMP = 20 # degC
    TIN_MIN_TEMP = 300 # degC, melting point = 232C

    base_energy = EP[:eTotalCapEnergy][end] * min_temp / (max_temp - min_temp)
    amb_energy = EP[:eTotalCapEnergy][end] * AMBIENT_TEMP / (max_temp - min_temp)
    min_energy = EP[:eTotalCapEnergy][end] * TIN_MIN_TEMP / (max_temp - min_temp)

    adjusted_self_disch = dfGen[end,:Self_Disch] / (1.0 + (min_temp - AMBIENT_TEMP) / (max_temp - min_temp))
    
    # Make sure we don't go to negative energy
    @constraint(EP, TEGS_melt_energy[t=1:myinputs["T"]], EP[:vS][end,t] >= min_energy - base_energy)

    @constraint(EP, tegs_cSoCBalStart[t in START_SUBPERIODS, tegs_rid], 
        EP[:vS][tegs_rid,t] 
        ==
		+ EP[:vS][tegs_rid,t+hours_per_subperiod-1]
        - (EP[:vP][tegs_rid,t] / dfGen[tegs_rid,:Eff_Down])
		+ (dfGen[tegs_rid,:Eff_Up] * EP[:vCHARGE][tegs_rid,t])
        - (adjusted_self_disch * (base_energy - amb_energy + EP[:vS][tegs_rid,t+hours_per_subperiod-1]))
    )

    @constraint(EP, tegs_cSoCBalInterior[t in INTERIOR_SUBPERIODS, tegs_rid], 
        EP[:vS][tegs_rid,t] 
        ==
		+ EP[:vS][tegs_rid,t-1]
        - (EP[:vP][tegs_rid,t] / dfGen[tegs_rid,:Eff_Down])
        + (dfGen[tegs_rid,:Eff_Up] * EP[:vCHARGE][tegs_rid,t])
        - (adjusted_self_disch * (base_energy - amb_energy + EP[:vS][tegs_rid,t-1]))
    )

    ## Constrain discharge power based on state of charge. Reaches rated capacity at 100% state of charge
    # {Discharge power} / {discharge capacity} <= {energy} / {energy capacity}
    # {Discharge power} / {discharge capacity} <= {energy} / {discharge capacity * dsratio}
    # {Discharge power} * {dsratio} <= {energy}
    # Regular constraint. Gives the correct energy balance, but vS can't be negative as vP >= 0
    @constraint(EP, TEGSeffpower[t=1:myinputs["T"]], EP[:vP][end,t] * dsratio <= EP[:vS][end,t])

    ## Constrain charging power based on state of charge. Reaches rated capacity at 0% state of charge
    # {Charge power} / {charge capacity} <= 1 - {energy} / {energy capacity}
    # {Charge power} / {charge capacity} <= 1 - {energy} / {charge capacity * csratio}
    # {Charge power} * {csratio} <= {charge capacity * csratio} - {energy}
    @constraint(EP, TEGSeffcharge[t=1:myinputs["T"]], EP[:vCHARGE][end,t] * csratio <= EP[:eTotalCapCharge][end] * csratio - EP[:vS][end,t])

    println("Solving Model")
    EP, solve_time = solve_model(EP, mysetup)
    myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

    # Run MGA if the MGA flag is set to 1 else only save the least cost solution
    println("Writing Output")
    elapsed_time = @elapsed write_outputs(EP, outputs_path, mysetup, myinputs)
    println("Time elapsed for writing is")
    println(elapsed_time)

    return objective_value(EP)

end

############################################
# Some necessary helper functions
############################################
include("../RunTools.jl")

function costpermwdischarge(T)
    T = Int(T)
    discharge_vs_temp = Dict{Float64, Float64}(
        2400 => 0.3097 * 1e6,
        2300 => 0.3427 * 1e6,
        2100 => 0.4926 * 1e6,
        1900 => 0.8924 * 1e6,
    )
    return discharge_vs_temp[T]
end

function costpermwcharge(T)
    return 0.0185 * 1e6
end

function costpermwhstor(T, lossrate)
    storage_vs_temp = Dict{Float64, Array{Float64}}(
        2400 => 1e3 * [15.827, 12.0636, 10.8529, 10.2619, 9.9158, 9.6914, 9.5363, 9.4243, 9.3411, 9.2779, 9.2292, 9.1915, 9.1622, 9.1395, 9.122],
        2300 => 1e3 * [15.1493, 11.7375, 10.6382, 10.1019, 9.7883, 9.5854, 9.4456, 9.3451, 9.2707, 9.2146, 9.1717, 9.1388, 9.1136, 9.0943, 9.0799],
        2100 => 1e3 * [13.8087, 11.0887, 10.2103, 9.7827, 9.5337, 9.3737, 9.2644, 9.1867, 9.1301, 9.0881, 9.0568, 9.0336, 9.0165, 9.0042, 8.9958],
        1900 => 1e3 * [12.4875, 10.4445, 9.7845, 9.4646, 9.2799, 9.1626, 9.0836, 9.0287, 8.9897, 8.9618, 8.9421, 8.9284, 8.9194, 8.9141, 8.9117],
    )
    return storage_vs_temp[T][floor(Int, lossrate)]
end

function annuitize(cost::Float64, lifetime::Float64, discount_rate::Float64)
    return cost * discount_rate / (1 - (1 + discount_rate)^(-lifetime))
end

function TEGS_costs(stortype::Int, T::Float64, lossrate::Float64, lifetime::Float64, discount_rate::Float64)
    discharge_cost = annuitize(costpermwdischarge(T), lifetime, discount_rate)
    storage_cost = annuitize(costpermwhstor(T, lossrate), lifetime, discount_rate)
    charge_cost = annuitize(costpermwcharge(T), lifetime, discount_rate)
    if stortype == 1
        costs = Dict{String, Float64}(
        "Inv_Cost_per_MWyr" => discharge_cost + charge_cost,
        "Fixed_OM_Cost_per_MWyr" => 0.25 * (discharge_cost + charge_cost),
        "Inv_Cost_per_MWhyr" => storage_cost,
        "Fixed_OM_Cost_per_MWhyr" => 0.25 * storage_cost,
    )
    elseif stortype == 2
        costs = Dict{String, Float64}(
        "Inv_Cost_per_MWyr" => discharge_cost,
        "Fixed_OM_Cost_per_MWyr" => 0.25 * discharge_cost,
        "Inv_Cost_per_MWhyr" => storage_cost,
        "Fixed_OM_Cost_per_MWhyr" => 0.25 * storage_cost,
        "Inv_Cost_Charge_per_MWyr" => charge_cost,
        "Fixed_OM_Cost_Charge_per_MWyr" => 0.25 * charge_cost,
    )
    else
        error("Invalid storage type")
    end
    return costs
end

function selectresource(dfGen::DataFrame, resourcename::String)
    dfgen_search = groupby(dfGen, :Resource)
    return dfgen_search[(resourcename,)]
end
    
function setTEGScosts!(dfGen::DataFrame, stortype::Int, T::Float64, lossrate::Float64, lifetime::Float64, discount_rate::Float64)
    TEGS_input = selectresource(dfGen, "TEGS")
    tegs_costs = TEGS_costs(stortype, T, lossrate, lifetime, discount_rate)
    for (param, value) in tegs_costs
        TEGS_input[!,param] .= value
    end   
end

############################################
# Case Definitions 
# All cases intended to be run from the run-file directory
############################################
root_dir = dirname(dirname(dirname(@__FILE__))) # Should be ../TEGS_runs
run_name = "temp_lossrate_sweep_stor2_v6_testcase"
dropbox_path = "/Users/rmacd/Dropbox/1_Academics/Research/22-TEGS_modelling/TEGS GenX shared folder/GenX_runs"
# dropbox_path = "D:/Dropbox/1_Academics/Research/22-TEGS_modelling/TEGS GenX shared folder/GenX_runs"
# dropbox_path = "/media/rmacd/LargeHD/Dropbox/1_Academics/Research/22-TEGS_modelling/TEGS GenX shared folder/GenX_runs"

location_dir = Dict{String, String}(
    "newEngland" => joinpath(root_dir, "data", "newEngland_stor2_v2"),
    "texas" => joinpath(root_dir, "data", "texas_stor2"),
)

logging_notes = Array{String, 1}()

emissions_levels = Dict{String, Float64}(
    "0.01"=>0.01
    )

emiss_targets = Dict{String, Dict{String, Float64}}(
    "newEngland" => Dict{String, Float64}(),
    "texas" => Dict{String, Float64}(),
)

for (loc_name, loc_emiss) in emiss_targets
    outputs_path = joinpath(dropbox_path, "outputs", loc_name, "emissions_and_baseline_v2", "baseline")
    # outputs_path = joinpath(root_dir, "outputs", loc_name, "emissions_and_baseline", "baseline")
    emissions_df = DataFrame(CSV.File(joinpath(outputs_path, "emissions.csv")))
    baseline_emiss = emissions_df.Total[1]
    for (emiss_name, emiss_level) in emissions_levels
        loc_emiss[emiss_name] = baseline_emiss * emiss_level
    end
end
        
temperatures = Array{Float64}([2400])
lossrates = Array{Float64}([3])

for (loc_name, loc_path) in location_dir
    case = loc_path
    genx_settings = get_settings_path(case, "genx_settings.yml") #Settings YAML file path
    mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters

    mysetup["CO2Cap"] = 1
    push!(logging_notes, "Emissions cap activated\n")

    scale_factor = mysetup["ParameterScale"] == 1 ? ModelScalingFactor : 1

    inputs_path = case
    settings_path = get_settings_path(case)

    # outputs_path_root = joinpath(root_dir, "outputs", loc_name, run_name)
    # outputs_path_root = joinpath(dropbox_path, "outputs", loc_name, run_name)
    outputs_path_root = dirname(@__FILE__)
    if !isdir(outputs_path_root)
        mkpath(outputs_path_root)
    end

    ### Cluster time series inputs if necessary and if specified by the user
    TDRpath = joinpath(case, mysetup["TimeDomainReductionFolder"])

    if mysetup["TimeDomainReduction"] == 1
        if !time_domain_reduced_files_exist(TDRpath)
            println("Clustering Time Series Data (Grouped)...")
            cluster_inputs(inputs_path, settings_path, mysetup)
        else
            println("Time Series Data Already Clustered.")
        end
    end

    ### Configure solver
    println("Configuring Solver")
    OPTIMIZER = configure_solver(mysetup["Solver"], settings_path)

    ### Load inputs
    println("Loading Inputs")
    myinputs = load_inputs(mysetup, case)

    LIFETIME = 30.
    DISCOUNT_RATE = 0.05
    STOR_TYPE = 2

    push!(logging_notes, "LIFETIME = $LIFETIME\n")
    push!(logging_notes, "DISCOUNT_RATE = $DISCOUNT_RATE\n")
    push!(logging_notes, "STOR_TYPE = $STOR_TYPE\n")
    
    #### Running a case
    for (emiss_name, emiss_target) in emiss_targets[loc_name]
        myinputs["dfMaxCO2"] = emiss_target * 1e3 / scale_factor
        for T in temperatures
            for lossrate in lossrates
                push!(logging_notes, "Running $(emiss_name)_$(T)_$(lossrate)_stor2 case\n")
                setTEGScosts!(myinputs["dfGen"], STOR_TYPE, T, lossrate, LIFETIME, DISCOUNT_RATE)
                TEGS_input = selectresource(myinputs["dfGen"], "TEGS")
                TEGS_input[!, "Self_Disch"] .= lossrate / 100 / 24. # Convert daily loss rate to hourly
                TEGS_input[!, "STOR"] .= 2

                # Calculate and save baseline emissions
                case_key = "$(emiss_name)_$(T)_$(lossrate)_stor2"
                dirlist = readdir(outputs_path_root)
                # exist_flag = false
                # for dir in dirlist
                #     if startswith(dir, string(case_key, "_ds"))
                #         println("$(dir) already exists. Skipping...")
                #         exist_flag = true
                #     end
                # end
                # if exist_flag
                #     continue
                # end
                outputs_path = joinpath(outputs_path_root, case_key)

                dsratio = 1.5
                csratio = 1.5
                max_temp = T
                min_temp = T - 500.0
                y = tegs_case(dsratio, csratio, max_temp, min_temp, 13, outputs_path, mysetup, myinputs, OPTIMIZER)

            end
        end
    end
end

# Write logging notes
open(joinpath(dirname(@__FILE__), "$(run_name).log"), "w") do f
    for line in logging_notes
        write(f, line)
    end
end