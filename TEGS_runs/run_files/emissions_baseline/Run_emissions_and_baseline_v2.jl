using GenX
using CSV
using DataFrames

############################################
# Some necessary helper functions
############################################
include("../RunTools.jl")

############################################
# Case Definitions 
# All cases intended to be run from the run-file directory
############################################
root_dir = dirname(dirname(dirname(@__FILE__))) # Should be ../TEGS_runs

location_dir = Dict{String, String}(
    "newEngland" => joinpath(root_dir, "data", "newEngland_stor2_v2"),
    # "texas" => joinpath(root_dir, "data", "texas"),
)

logging_notes = Array{String, 1}()

for (loc_name, loc_path) in location_dir
    run_name = "emissions_and_baseline_v2"
    dropbox_path = "D:/Dropbox/1_Academics/Research/22-TEGS_modelling/TEGS GenX shared folder/GenX_runs"

    push!(logging_notes, "Running $loc_name case\n")
    case = loc_path
    genx_settings = get_settings_path(case, "genx_settings.yml") #Settings YAML file path
    mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters

    inputs_path = case
    settings_path = get_settings_path(case)

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
    println("Configuring Solver")
    OPTIMIZER = configure_solver(mysetup["Solver"], settings_path)

    #### Running a case

    ### Load inputs
    println("Loading Inputs")
    myinputs = load_inputs(mysetup, case)

    # Calculate and save baseline emissions
    outputs_path = joinpath(outputs_path_root, "baseline")
    EP = build_solve_write(outputs_path, mysetup, myinputs, OPTIMIZER)

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

# Write logging notes
open(joinpath(dirname(@__FILE__), "emissions_and_baseline.log"), "w") do f
    for line in logging_notes
        write(f, line)
    end
end