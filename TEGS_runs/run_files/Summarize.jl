using CSV
using DataFrames

############################################

root_dir = dirname(dirname(@__FILE__)) # Should be ../TEGS_runs
outputs_dir = joinpath(root_dir, "outputs")

location_dir = Dict{String, String}(
    "newEngland" => joinpath(root_dir, "data", "newEngland"),
    # "texas" => joinpath(root_dir, "data", "texas"),
)
cases = ["emissions_and_baselineTest"]

# Info to summarize for each resource
resource_cols = Dict{String, Array{Symbol}}(
    "capacity.csv" => [
        :Resource,
        :EndCap, 
        :EndEnergyCap, 
        :EndChargeCap
    ],
)

zone_cols = Dict{String, Array{Symbol}}()

files_to_search = union(collect(keys(resource_cols)), collect(keys(zone_cols)))

# resource_summ = DataFrame()
resource_summ = Dict{String, Dict{String, Any}}()

for (loc_name, loc_path) in location_dir
    for case in cases
        result_dir = joinpath(outputs_dir, loc_name, case)
        println("Printing results from $result_dir")

        results = readdir(result_dir)
        for result in results
            for result_file in files_to_search
                capacity_df = DataFrame(CSV.File(joinpath(result_dir, result, result_file)))
                EndCap = Dict{String, Array{Float64}}(
                    "EndCap" => Dict{String, Array{Float64}}(capacity_df.Resource .=> capacity_df.EndCap),
                    "EndEnergyCap" => capacity_df.EndEnergyCap,
                    "EndChargeCap" => capacity_df.EndChargeCap,
                )
                
                # summ_cols = intersect(propertynames(capacity_df), resource_cols[result_file])
                # println(summ_cols)
                # summ_cols = setdiff(summ_cols, propertynames(resource_summ))
                # println(summ_cols)
                # resource_summ = hcat(resource_summ, capacity_df[!, summ_cols], makeunique=true)
            end
        end
    end
end

