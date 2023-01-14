using CSV
using JSON3
using DataFrames

include("Summarize.jl")

# Find all summary files below this
root_dir = "/Users/rmacd/Dropbox/1_Academics/Research/22-TEGS_modelling/TEGS GenX shared folder/GenX_runs"
outputs_dir = joinpath(root_dir, "outputs")

location_dir = Dict{String, String}(
    "newEngland" => joinpath(outputs_dir, "newEngland"),
    "texas" => joinpath(outputs_dir, "texas"),
)

summ_paths = Dict{String, String}()
for (loc_name, loc_path) in location_dir
    getlocsummaryfiles!(summ_paths, loc_path, loc_name)
end

json_string = read(summ_paths["newEngland_emissions_and_baseline"], String)
example_result = JSON3.read(json_string, Dict)