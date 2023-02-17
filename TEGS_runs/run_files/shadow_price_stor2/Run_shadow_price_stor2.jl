using GenX
using CSV
using DataFrames
using JuMP

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

function setFreeTEGS(dfGen::DataFrame, stortype::Int, T::Float64, lossrate::Float64, lifetime::Float64, discount_rate::Float64)
    if stortype == 1
        costs = Dict{String, Float64}(
        "Inv_Cost_per_MWyr" => 0,
        "Fixed_OM_Cost_per_MWyr" => 0,
        "Inv_Cost_per_MWhyr" => 0,
        "Fixed_OM_Cost_per_MWhyr" => 0,
    )
    elseif stortype == 2
        costs = Dict{String, Float64}(
        "Inv_Cost_per_MWyr" => 0,
        "Fixed_OM_Cost_per_MWyr" => 0,
        "Inv_Cost_per_MWhyr" => 0,
        "Fixed_OM_Cost_per_MWhyr" => 0,
        "Inv_Cost_Charge_per_MWyr" => 0,
        "Fixed_OM_Cost_Charge_per_MWyr" => 0,
    )
    end
    TEGS_input = selectresource(dfGen, "TEGS")
    costs = TEGS_costs(stortype, T, lossrate, lifetime, discount_rate)
    for (param, value) in costs
        TEGS_input[!,param] .= value
    end
end

############################################
# Case Definitions 
# All cases intended to be run from the run-file directory
############################################
root_dir = dirname(dirname(dirname(@__FILE__))) # Should be ../TEGS_runs
run_name = "shadow_price_stor2"
dropbox_path = "D:/Dropbox/1_Academics/Research/22-TEGS_modelling/TEGS GenX shared folder/GenX_runs"
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
    emissions_df = DataFrame(CSV.File(joinpath(outputs_path, "emissions.csv")))
    baseline_emiss = emissions_df.Total[1]
    for (emiss_name, emiss_level) in emissions_levels
        loc_emiss[emiss_name] = baseline_emiss * emiss_level
    end
end
        
# temperatures = Array{Float64}([2400, 2300, 2100, 1900])
# lossrates = Array{Float64}([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
temperatures = Array{Float64}([2400])
lossrates = Array{Float64}([3])

TEGSdischargecaps = Array{Float64}([10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000])
TEGSdischarge_chargeratio = Array{Float64}([1, 2, 3, 4, 5, 6])
TEGSdischarge_storageratio = Array{Float64}([10, 20, 30, 40, 50, 60])

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
    outputs_path_root = joinpath(dropbox_path, "outputs", loc_name, run_name)
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
    OPTIMIZER = configure_solver(mysetup["Solver"], settings_path)

    ### Load inputs
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
                for d_c_ratio in TEGSdischarge_chargeratio
                    for d_s_ratio in TEGSdischarge_storageratio
                        for discap in TEGSdischargecaps
                            sub_case_name = "$(emiss_name)_$(T)_$(lossrate)_$(discap)_$(d_c_ratio)_$(d_s_ratio)_stor2"
                            outputs_path = joinpath(outputs_path_root, sub_case_name)
                            if isdir(outputs_path)
                                push!(logging_notes, "Skipping $(sub_case_name) case -- already exists\n")
                                continue
                            end
                            push!(logging_notes, "Running $(sub_case_name) case\n")
                            setFreeTEGS(myinputs["dfGen"], STOR_TYPE, T, lossrate, LIFETIME, DISCOUNT_RATE)
                            TEGS_input = selectresource(myinputs["dfGen"], "TEGS")
                            TEGS_input[!, "Self_Disch"] .= lossrate / 100 / 24. # Convert daily loss rate to hourly
                            TEGS_input[!, "STOR"] .= 2

                            # Calculate and save baseline emissions
                            EP = generate_model(mysetup, myinputs, OPTIMIZER)

                            @constraint(EP, tegsdiscap, EP[:eTotalCap][end] == discap)
                            @constraint(EP, tegscharcap, EP[:eTotalCapCharge][end] == discap * d_c_ratio)
                            @constraint(EP, tegsstorcap, EP[:eTotalCapEnergy][end] == discap * d_s_ratio)

                            EP, solve_time = solve_model(EP, mysetup)
                            myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

                            # Run MGA if the MGA flag is set to 1 else only save the least cost solution
                            elapsed_time = @elapsed write_outputs(EP, outputs_path, mysetup, myinputs)

                            open(joinpath(outputs_path, "shadow_price_summ.csv"), "w") do f
                                println(f, "name, discharge, charge, storage")
                                println(f, "target cap, $(discap), $(discap * d_c_ratio), $(discap * d_s_ratio)")
                                println(f, "end cap, $(value(EP[:eTotalCap][end])), $(value(EP[:eTotalCapCharge][end])), $(value(EP[:eTotalCapEnergy][end]))")
                                println(f, "value, $(dual.(EP[:tegsdiscap])), $(dual.(EP[:tegscharcap])), $(dual.(EP[:tegsstorcap]))")
                            end
                            
                            push!(logging_notes, "$(sub_case_name) results written to $outputs_path\n")

                            if any([iszero.(dual.(EP[:tegsdiscap])), iszero.(dual.(EP[:tegscharcap])), iszero.(dual.(EP[:tegsstorcap]))])
                                break
                            end
                        end
                    end
                end
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