using CSV
using DataFrames

############################################
# Helper functions
############################################
function numresources(filename::String)
    """Return the number of resources in a file"""
    csvfile = CSV.File(filename)
    if :Resource in propertynames(csvfile)
        return length(csvfile.Resource)
    else
        error("$filename does not have a Resource column") 
    end
end

function checknumresources(result_dir::String, result_file::String="capacity.csv")
    """Return the number of resources in a system, based off the specified file"""
    result = readdir(result_dir)[1]
    testfilename = joinpath(result_dir, result, result_file)
    if isfile(testfilename)
        return numresources(testfilename)
    else
        error("Could not check the number of resources\nThis file does not exist: $testfilename")
    end
end



############################################
# Main Loop
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
        "Resource",
        "EndCap", 
        "EndEnergyCap", 
        "EndChargeCap"
    ],
)

zone_cols = Dict{String, Array{Symbol}}()

files_to_search = union(collect(keys(resource_cols)), collect(keys(zone_cols)))

# resource_summ = DataFrame()
resource_summ = Dict{String, Dict{String, Any}}()

# We'll use dictionaries to store everything to reduce the chance of
# mislabeling data due to differences in the order of the columns
# We'll take a slight performance hit, but we don't have much data
# This will also make it easier to write to a JSON file

# For parameter / output of interest, e.g. end capacity, annual generation, etc.
# We'll create a dictionary of the different resources,
# and each entry will be a dictionary of the different cases

# Summary stats of the whole system will be in a separate dictionary

resourec_summ = Dict{String, Any}()

for (loc_name, loc_path) in location_dir
    for case in cases
        result_dir = joinpath(outputs_dir, loc_name, case)
        println("Printing results from $result_dir")

        results = readdir(result_dir)
        for result in results
            for result_file in files_to_search
                resultdata_df = DataFrame(CSV.File(joinpath(result_dir, result, result_file)))
                summ_cols = intersect(names(resultdata_df), resource_cols[result_file])
                summ_cols = setdiff(summ_cols, collect(keys(resource_summ)))
                for param in summ_cols
                    for resource in resultdata_df.Resource
                        resource_summ[param][resource] = Dict{String, Float64}(results => )
                    end
                    resource_summ[param] = Dict{String, Any}()
                end

                EndCap = Dict{String, Any}(
                    "EndCap" => Dict{String, Dict{String, Float64}}(capacity_df.Resource .=> Dict{String, Float64}(result => capacity_df.EndCap)),
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

