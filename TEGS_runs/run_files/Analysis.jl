using CSV
using JSON3
using DataFrames

function resourcesumm2df(resource_summ::Dict{String, Any})
    """Convert a resource summary to a dataframe"""
    df_inputs =Dict{String, Any}(collect(keys(resource_summ)) .=> collect.(values.(collect(values(resource_summ)))))
    df_inputs["Resource"] = collect(keys(collect(values(resource_summ))[1]))
    return select!(DataFrame(df_inputs), :Resource, Not(:Resource)) # Reorder the columns so Resource is first
end

function resourcesumm2df(resource_summ_filename::String)
    """Convert a resource summary JSON file to a dataframe"""
    summ_json_string = read(resource_summ_filename, String)
    resource_summ = JSON3.read(summ_json_string, Dict)
    return resourcesumm2df(resource_summ)
end

# Test path from "...TEGS_runs/run_files"
summ_path = "../outputs/newEngland/emissions_and_baseline/emissions_and_baseline_summary.json"

resource_summ_df = resourcesumm2df(summ_path)