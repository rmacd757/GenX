using GenX
using CSV
using DataFrames

############################################
# Some necessary helper functions
############################################
include("../RunTools.jl")

function costpermwdischarge(T)
    discharge_vs_temp = Dict{Float64, Float64}(
        2400 => 0.3097,
        2300 => 0.3427,
        2100 => 0.4926,
        1900 => 0.8924,
    )
    return discharge_vs_temp[T]
end

function costpermwcharge(T)
    return 0.0185

function costpermwhstor(T, lossrate)
    storage_vs_temp = Dict{Float64, Array{Float64}}(
        2400 => [15.827, 12.0636, 10.8529, 10.2619, 9.9158, 9.6914, 9.5363, 9.4243, 9.3411, 9.2779, 9.2292, 9.1915, 9.1622, 9.1395, 9.122],
        2300 => [15.1493, 11.7375, 10.6382, 10.1019, 9.7883, 9.5854, 9.4456, 9.3451, 9.2707, 9.2146, 9.1717, 9.1388, 9.1136, 9.0943, 9.0799],
        2100 => [13.8087, 11.0887, 10.2103, 9.7827, 9.5337, 9.3737, 9.2644, 9.1867, 9.1301, 9.0881, 9.0568, 9.0336, 9.0165, 9.0042, 8.9958],
        1900 => [12.4875, 10.4445, 9.7845, 9.4646, 9.2799, 9.1626, 9.0836, 9.0287, 8.9897, 8.9618, 8.9421, 8.9284, 8.9194, 8.9141, 8.9117],
    )
    return storage_vs_temp[T][lossrate]
end

############################################
# Case Definitions 
# All cases intended to be run from the run-file directory
############################################
root_dir = dirname(dirname(dirname(@__FILE__))) # Should be ../TEGS_runs

location_dir = Dict{String, String}(
    "newEngland" => joinpath(root_dir, "data", "newEngland"),
    "texas" => joinpath(root_dir, "data", "texas"),
)

logging_notes = Array{String, 1}()

for (loc_name, loc_path) in location_dir
    case = loc_path
    genx_settings = get_settings_path(case, "genx_settings.yml") #Settings YAML file path
    mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters

    inputs_path = case
    settings_path = get_settings_path(case)

    outputs_path_root = joinpath(root_dir, "outputs", loc_name, "emissions_and_baseline")
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

    #### Running a case

    ### Load inputs
    println("Loading Inputs")
    myinputs = load_inputs(mysetup, case)

    # Calculate and save baseline emissions
    outputs_path = joinpath(outputs_path_root, "baseline")
    EP = build_solve_write(outputs_path, mysetup, myinputs, OPTIMIZER)

    println
    push!(logging_notes, "Baseline results written to $outputs_path\n")

    # Find the emissions level
    # Load from the written results to make this usable with old results
    # DataFrame is overkill here, but more flexible
    emissions_df = DataFrame(CSV.File(joinpath(outputs_path, "emissions.csv")))
    baseline_emiss = emissions_df.Total[1]

    push!(logging_notes, "Baseline emissions were $baseline_emiss\n")

    emissions_levels = Dict{String, Float64}(
        "0.20"=>0.2,
        "0.10"=>0.1, 
        "0.05"=>0.05,
        "0.01"=>0.01)
    scale_factor = mysetup["ParameterScale"] == 1 ? ModelScalingFactor : 1

    
    mysetup["CO2Cap"] = 1
    myinputs = load_inputs(mysetup, case) # Reload CO2_cap settings
    push!(logging_notes, "Emissions cap activated and inputs reloaded\n")

    for (emiss_name, emiss_level) in emissions_levels
        emiss_target = emiss_level * baseline_emiss
        push!(logging_notes, "Emissions target is $emiss_target\n")

        myinputs["dfMaxCO2"] = emiss_target * 1e3 / scale_factor

        # Calculate and save baseline emissions
        outputs_path = joinpath(outputs_path_root, emiss_name)
        EP = build_solve_write(outputs_path, mysetup, myinputs, OPTIMIZER)

        push!(logging_notes, "$emiss_target results written to $outputs_path\n")
    end
end