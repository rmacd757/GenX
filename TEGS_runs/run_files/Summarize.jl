using CSV
using JSON3
using DataFrames

############################################
# Helper functions
############################################
function shortdirnames(dir::String)
    """Return the names of all subdirectories in a directory"""
    return basename.(filter(isdir, readdir(dir, join=true)))
end

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
    result = shortdirnames(result_dir)[1]
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

function summarizefile!(
    resource_summ::Dict{String, Any},
    result_filepath::String,
    resource_cols::Dict{String, Array{String}},
    zone_cols::Dict{String, Array{String}},
    )
    """
    Summarize the results of a single file\n
    Summarizes the parameters in resource_cols and zone_cols
    from the file at result_filepath into resource_summ\n
    """
    result_data = CSV.File(result_filepath)
    path_comp = splitpath(result_filepath)
    result_name = path_comp[end-1]
    result_file = path_comp[end] # aka: basename(result_filepath)
    summ_params = intersect(string.(propertynames(result_data)), resource_cols[result_file])
    for param in summ_params
        for (idx, resource) in enumerate(result_data.Resource)
            resource_summ[param][resource][result_name] = result_data[param][idx]
        end
    end
end

function summarizerun(
    result_dir::String, 
    resource_cols::Dict{String, Array{String}},
    zone_cols::Dict{String, Array{String}},
    )
    """Summarize all the results from the series of optimizations in result_dir"""
    println("Printing results from $result_dir")

    files_to_search = union(collect(keys(resource_cols)), collect(keys(zone_cols)))

    param_names = String[]
    for params in values(resource_cols)
        param_names = union(param_names, params)
    end
    resource_summ = initparamresults(param_names)

    resource_names = getresourcenames(result_dir)
    for param in keys(resource_summ)
        resource_summ[param] = initresourceresults(resource_names)
    end

    case_results = shortdirnames(result_dir)
    for result in case_results
        for result_file in files_to_search
            summarizefile!(
                resource_summ, 
                joinpath(result_dir, result, result_file), 
                resource_cols, 
                zone_cols)
        end
    end
    return resource_summ
end

function summarizerunandsave(
    result_dir::String, 
    resource_cols::Dict{String, Array{String}}, 
    zone_cols::Dict{String, Array{String}}
    )
    """Summarize all the results from the series of optimizations in result_dir and save them to a JSON file"""
    resource_summ = summarizerun(result_dir, resource_cols, zone_cols)
    case_name = basename(result_dir)
    open(joinpath(result_dir, "$(case_name)_summary.json"), "w") do io
        JSON3.pretty(io, resource_summ)
    end
end

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

############################################
# Example Main Loop
############################################
# root_dir = dirname(dirname(@__FILE__)) # Should be ../TEGS_runs
# outputs_dir = joinpath(root_dir, "outputs")

# location_dir = Dict{String, String}(
#     "newEngland" => joinpath(root_dir, "data", "newEngland"),
#     # "texas" => joinpath(root_dir, "data", "texas"),
# )
# cases = ["emissions_and_baseline"]

# # Info to summarize for each resource
# resource_cols = Dict{String, Array{String}}(
#     "capacity.csv" => [
#         "EndCap", 
#         "EndEnergyCap", 
#         "EndChargeCap"
#     ],
# )

# zone_cols = Dict{String, Array{String}}()

# for (loc_name, loc_path) in location_dir
#     for case in cases
#         result_dir = joinpath(outputs_dir, loc_name, case)
#         summarizerunandsave(result_dir, resource_cols, zone_cols)
#     end
# end




