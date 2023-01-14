using CSV
using JSON3
using DataFrames

include("../Summarize.jl")

# root_dir = dirname(dirname(dirname(@__FILE__))) # Should be ../TEGS_runs
root_dir = "/Users/rmacd/Dropbox/1_Academics/Research/22-TEGS_modelling/TEGS GenX shared folder/GenX_runs"
outputs_dir = joinpath(root_dir, "outputs")

location_dir = Dict{String, String}(
    "newEngland" => joinpath(root_dir, "outputs", "newEngland"),
    "texas" => joinpath(root_dir, "outputs", "texas"),
)
cases = ["emissions_and_baseline", "temp_lossrate_sweep", "temp_lossrate_sweep_stor2"]

# Info to summarize for each resource
resource_cols = Dict{String, Array{String}}(
    "capacity.csv" => [
        "EndCap", 
        "EndEnergyCap", 
        "EndChargeCap"
    ],
)

zone_cols = Dict{String, Array{String}}(
    "costs.csv" => [
        "cTotal"
    ]
)

summ_paths = Dict{String, String}()

for (loc_name, loc_path) in location_dir
    for case in cases
        result_dir = joinpath(outputs_dir, loc_name, case)
        summ_file_name = summarizerunandsave(result_dir, resource_cols, zone_cols)
        summ_paths[string(loc_name, "_", case)] = summ_file_name
    end
end

# Test path from "...TEGS_runs/run_files"
# summ_paths = Dict{String, String}(
#     "newEngland" => "../../outputs/newEngland/emissions_and_baseline/emissions_and_baseline_summary.json",
#     "texas" => "../../outputs/texas/emissions_and_baseline/emissions_and_baseline_summary.json",
# )

# resource_summ_df = resourcesumm2df(summ_path)