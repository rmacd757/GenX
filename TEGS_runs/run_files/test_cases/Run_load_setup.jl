using GenX
using DataFrames

############################################
# Some necessary helper functions
############################################

function get_settings_path(case::AbstractString)
    return joinpath(case, "Settings")
end

function get_settings_path(case::AbstractString, filename::AbstractString)
    return joinpath(get_settings_path(case), filename)
end

function get_default_output_folder(case::AbstractString)
    return joinpath(case, "Results")
end

function time_domain_reduced_files_exist(tdrpath)
    tdr_load = isfile(joinpath(tdrpath,"Load_data.csv"))
    tdr_genvar = isfile(joinpath(tdrpath,"Generators_variability.csv"))
    tdr_fuels = isfile(joinpath(tdrpath,"Fuels_data.csv"))
    return (tdr_load && tdr_genvar && tdr_fuels)
end

############################################
# Case Definitions 
# All cases intended to be run from the run-file directory
############################################
root_dir = dirname(dirname(dirname(@__FILE__))) # Should be ../TEGS_runs

location_dir = Dict{String, String}(
    "newEngland" => joinpath(root_dir, "data", "newEngland"),
    # "texas" => joinpath(root_dir, "data", "texas"),
)
loc_name = "newEngland"
loc_path = joinpath(root_dir, "data", "newEngland")

scenarios = ["results"]
scenario = "results"

case = loc_path
genx_settings = get_settings_path(case, "genx_settings.yml") #Settings YAML file path
mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters

inputs_path = case
settings_path = get_settings_path(case)

outputs_path = joinpath(root_dir, "outputs", "tests", loc_name)
outputs_path = joinpath(outputs_path, scenario)

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

# println("Generating the Optimization Model")
# EP = generate_model(mysetup, myinputs, OPTIMIZER)

# println("Solving Model")
# EP, solve_time = solve_model(EP, mysetup)
# myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

# # Run MGA if the MGA flag is set to 1 else only save the least cost solution
# println("Writing Output")
# elapsed_time = @elapsed write_outputs(EP, outputs_path, mysetup, myinputs)
# println("Time elapsed for writing is")
# println(elapsed_time)