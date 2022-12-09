using GenX
using DataFrames

############################################

# Some necessary helper functions
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

case = dirname(@__FILE__)
genx_settings = get_settings_path(case, "genx_settings.yml") #Settings YAML file path
mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters

inputs_path = case
settings_path = get_settings_path(case)

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

invcostarr = (6000.0:1000.0:20000.0) # This should be the annuitized cost [$/MWe], NOT the overnight cost

for idx in 1:length(invcostarr)

    # This changes the investment cost of the last generator in your Generators_data.csv file
    # Change the column name if you want to change a different file
    # Change the index if you want to change a different generator, e.g. [end-1] = the second last
    # This is not a robust solution, but should work for your purposes
    myinputs["dfGen"][!,"Inv_cost_per_MWyr"][end] = invcostarr[idx]

    println("Generating the Optimization Model")
    EP = generate_model(mysetup, myinputs, OPTIMIZER)

    println("Solving Model")
    EP, solve_time = solve_model(EP, mysetup)
    myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

    # Run MGA if the MGA flag is set to 1 else only save the least cost solution
    println("Writing Output")
    outputs_path = get_default_output_folder(case)
    elapsed_time = @elapsed write_outputs(EP, outputs_path, mysetup, myinputs)
    println("Time elapsed for writing is")
    println(elapsed_time)
    if mysetup["ModelingToGenerateAlternatives"] == 1
        println("Starting Model to Generate Alternatives (MGA) Iterations")
        mga(EP, case, mysetup, myinputs, outputs_path)
    end

    if mysetup["MethodofMorris"] == 1
        println("Starting Global sensitivity analysis with Method of Morris")
        morris(EP, inputs_path, mysetup, myinputs, outputs_path, OPTIMIZER)
    end
end