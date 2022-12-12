using CSV
using JSON3
using DataFrames

############################################
# Helper functions
############################################
function resourcenames(filename::String)
    """Return the number of resources in a file"""
    csvfile = CSV.File(filename)
    if :Resource in propertynames(csvfile)
        return csvfile.Resource
    else
        error("$filename does not have a Resource column") 
    end
end

function getresourcenames(result_dir::String, result_file::String="capacity.csv")
    """Return the number of resources in a system, based off the specified file"""
    result = readdir(result_dir)[1]
    testfilename = joinpath(result_dir, result, result_file)
    if isfile(testfilename)
        return resourcenames(testfilename)
    else
        error("Could not check the number of resources\nThis file does not exist: $testfilename")
    end
end

function initresourceresults(resource_names::Array{String})
    """Create a dictionary of the results for each resource"""
    resource_results = Dict{String, Any}()
    for resource in resource_names
        resource_results[resource] = Dict{String, Float64}()
    end
    return resource_results
end

function initparamresults(param_names::Array{String})
    """Create a dictionary of the results for each parameter"""
    param_results = Dict{String, Any}()
    for param in param_names
        param_results[param] = Dict{String, Any}()
    end
    return param_results
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
cases = ["emissions_and_baseline"]

# Info to summarize for each resource
resource_cols = Dict{String, Array{String}}(
    "capacity.csv" => [
        "EndCap", 
        "EndEnergyCap", 
        "EndChargeCap"
    ],
)

zone_cols = Dict{String, Array{Symbol}}()

files_to_search = union(collect(keys(resource_cols)), collect(keys(zone_cols)))

# resource_summ = DataFrame()
param_names = String[]
for params in values(resource_cols)
    global param_names = union(param_names, params)
end
resource_summ = initparamresults(param_names)

# We'll use dictionaries to store everything to reduce the chance of
# mislabeling data due to differences in the order of the columns
# We'll take a slight performance hit, but we don't have much data
# This will also make it easier to write to a JSON file

# For parameter / output of interest, e.g. end capacity, annual generation, etc.
# We'll create a dictionary of the different resources,
# and each entry will be a dictionary of the different cases

# Summary stats of the whole system will be in a separate dictionary

for (loc_name, loc_path) in location_dir
    for case in cases
        result_dir = joinpath(outputs_dir, loc_name, case)
        println("Printing results from $result_dir")

        resource_names = getresourcenames(result_dir)
        for param in keys(resource_summ)
            resource_summ[param] = initresourceresults(resource_names)
        end

        case_results = basename.(filter(isdir, readdir(result_dir, join=true)))
        for result in case_results
            for result_file in files_to_search
                result_data = CSV.File(joinpath(result_dir, result, result_file))

                summ_params = intersect(string.(propertynames(result_data)), resource_cols[result_file])

                for param in summ_params
                    for (idx, resource) in enumerate(result_data.Resource)
                        resource_summ[param][resource][result] = result_data[param][idx]
                    end
                end
            end
        end
        open(joinpath(result_dir, "$(case)_summary.json"), "w") do io
            JSON3.pretty(io, resource_summ)
        end
    end
end



