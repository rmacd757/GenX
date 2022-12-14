using CSV
using JSON3
using DataFrames

include("../Summarize.jl")

oot_dir = dirname(dirname(@__FILE__)) # Should be ../TEGS_runs
outputs_dir = joinpath(root_dir, "outputs")

location_dir = Dict{String, String}(
    "newEngland" => joinpath(root_dir, "data", "newEngland"),
    "texas" => joinpath(root_dir, "data", "texas"),
)
cases = ["emissions_and_baseline"]

# Info to summarize for each resource
resource_cols = Dict{String, Array{String}}(
    "capacity.csv" => [
        "EndCap", 
        "EndEnergyCap", 
        "EndChargeCap"
    ],
)

zone_cols = Dict{String, Array{String}}()

for (loc_name, loc_path) in location_dir
    for case in cases
        result_dir = joinpath(outputs_dir, loc_name, case)
        summarizerunandsave(result_dir, resource_cols, zone_cols)
    end
end

# Test path from "...TEGS_runs/run_files"
summ_paths = Dict{String, String}(
    "newEngland" => "../../outputs/newEngland/emissions_and_baseline/emissions_and_baseline_summary.json",
    "texas" => "../../outputs/texas/emissions_and_baseline/emissions_and_baseline_summary.json",
)

# resource_summ_df = resourcesumm2df(summ_path)