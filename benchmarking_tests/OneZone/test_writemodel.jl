using GenX
using JuMP
using OrderedCollections
using DataFrames
using CSV
using BenchmarkTools

## Tell GenX where the inputs and ouputs are located
case_path = dirname(@__FILE__)
inputs_path = case_path

## Load helper functions to let us run the cases
include(joinpath(dirname(case_path),"run_helpers.jl"))

## Load settings
genx_settings = get_settings_path(inputs_path, "genx_settings.yml") #Settings YAML file path
mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters
settings_path = get_settings_path(inputs_path)

## Create and configure solver
println("Configuring Solver")
OPTIMIZER = configure_solver(mysetup["Solver"], settings_path)

#### Running a case

## Load inputs
println("Loading Inputs")
myinputs = load_inputs(mysetup, inputs_path)

outputs_path = joinpath(case_path, "Results")
mkpath(dirname(outputs_path))
    
## Generate model
println("Generating the Optimization Model")
EP = generate_model(mysetup, myinputs, OPTIMIZER)

## Solve model
println("Solving Model")
EP, solve_time = solve_model(EP, mysetup)
myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

## Run MGA if the MGA flag is set to 1 else only save the least cost solution
println("Writing Output")
# outputs_path = get_default_output_folder(outputs_path)

## Write outputs
benchmark_result = @benchmark write_outputs(EP, outputs_path, mysetup, myinputs) seconds=30.0 samples=50 evals=1

# Return the time (nanoseconds), garbage collection time (nanoseconds), memory (megabytes), and number of allocations
df = DataFrame(
    times_ns = benchmark_result.times, 
    gctimes_ns = benchmark_result.gctimes,
    memory_mb = benchmark_result.memory ./ 1e6,
    allocs = benchmark_result.allocs
)

# Write the results to a CSV file
if !isdir(outputs_path)
    mkpath(outputs_path)
end
CSV.write(joinpath(outputs_path, "writing_model_benchmark_results.csv"), df)