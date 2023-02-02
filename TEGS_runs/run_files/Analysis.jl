using CSV
using JSON3
using DataFrames
using Plots; pythonplot()

include("Summarize.jl")

# Find all summary files below this
dropbox_dir = "/Users/rmacd/Dropbox/"
# dropbox_dir = "D:/Dropbox"
root_dir = joinpath(dropbox_dir, "1_Academics/Research/22-TEGS_modelling/TEGS GenX shared folder/GenX_runs")
outputs_dir = joinpath(root_dir, "outputs")

location_dir = Dict{String, String}(
    "newEngland" => joinpath(outputs_dir, "newEngland"),
    # "texas" => joinpath(outputs_dir, "texas"),
)

summ_paths = Dict{String, String}()
for (loc_name, loc_path) in location_dir
    getlocsummaryfiles!(summ_paths, loc_path, loc_name)
end

for loc_name in collect(keys(location_dir))
    json_string = read(summ_paths["$(loc_name)_emissions_and_baseline_v2"], String)
    example_result = JSON3.read(json_string, Dict)

    json_string = read(summ_paths["$(loc_name)_temp_lossrate_sweep"], String)
    example_result1 = JSON3.read(json_string, Dict)

    json_string = read(summ_paths["$(loc_name)_temp_lossrate_sweep_stor2_v2"], String)
    example_result2 = JSON3.read(json_string, Dict)

    # function getmatch(result, idx)
    #     return split(result, "_")[idx]
    # end

    emiss_levels = Array{String}(["0.01", "0.05", "0.10", "0.20"])
    temperatures = Array{Float64}([2400, 2300, 2100, 1900])
    lossrates = Array{Float64}([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])

    plots = []

    for emiss_level in emiss_levels
        temp2 = Array{Float64}(undef, length(temperatures), length(lossrates))
        for (i, temp_val) in enumerate(temperatures)
            for (j, lossrate_val) in enumerate(lossrates)
                temp2[i,j] = example_result2["cTotal"]["Total"]["$(emiss_level)_$(temp_val)_$(lossrate_val)_stor2"]
            end
        end

        temp1 = Array{Float64}(undef, length(temperatures), length(lossrates))
        for (i, temp_val) in enumerate(temperatures)
            for (j, lossrate_val) in enumerate(lossrates)
                temp1[i,j] = example_result1["cTotal"]["Total"]["$(emiss_level)_$(temp_val)_$(lossrate_val)"]
            end
        end

        temp0 = Array{Float64}(undef, length(temperatures), length(lossrates))
        for (i, temp_val) in enumerate(temperatures)
            for (j, lossrate_val) in enumerate(lossrates)
                temp0[i,j] = example_result["cTotal"]["Total"]["$(emiss_level)"]
            end
        end

        cap2 = Array{Float64}(undef, length(temperatures), length(lossrates))
        for (i, temp_val) in enumerate(temperatures)
            for (j, lossrate_val) in enumerate(lossrates)
                cap2[i,j] = example_result2["EndCap"]["TEGS"]["$(emiss_level)_$(temp_val)_$(lossrate_val)_stor2"]
            end
        end

        cap1 = Array{Float64}(undef, length(temperatures), length(lossrates))
        for (i, temp_val) in enumerate(temperatures)
            for (j, lossrate_val) in enumerate(lossrates)
                cap1[i,j] = example_result1["EndCap"]["TEGS"]["$(emiss_level)_$(temp_val)_$(lossrate_val)"]
            end
        end

        cap0 = Array{Float64}(undef, length(temperatures), length(lossrates))
        for (i, temp_val) in enumerate(temperatures)
            for (j, lossrate_val) in enumerate(lossrates)
                cap0[i,j] = example_result["EndCap"]["TEGS"]["$(emiss_level)"]
            end
        end

        char2 = Array{Float64}(undef, length(temperatures), length(lossrates))
        for (i, temp_val) in enumerate(temperatures)
            for (j, lossrate_val) in enumerate(lossrates)
                char2[i,j] = example_result2["EndChargeCap"]["TEGS"]["$(emiss_level)_$(temp_val)_$(lossrate_val)_stor2"]
            end
        end

        char1 = Array{Float64}(undef, length(temperatures), length(lossrates))
        for (i, temp_val) in enumerate(temperatures)
            for (j, lossrate_val) in enumerate(lossrates)
                char1[i,j] = example_result1["EndCap"]["TEGS"]["$(emiss_level)_$(temp_val)_$(lossrate_val)"]
            end
        end

        c1 = contour(lossrates, temperatures, -100 .* [(temp1.-temp0)./temp0,(temp2.-temp0)./temp0], clabels=true, color=:turbo)
        c2 = contour(lossrates, temperatures, [cap1, cap2], clabels=true, color=:turbo, title="$(emiss_level) D")
        c3 = contour(lossrates, temperatures, [char2], clabels=true, color=:turbo, title="$(emiss_level) C")
        c4 = contour(lossrates, temperatures, [cap2./char2], clabels=true, color=:turbo, title="$(emiss_level) D / C")
        push!(plots, c1)
        push!(plots, c2)
        push!(plots, c3)
        push!(plots, c4)
    end
    plot(plots..., layout=(length(emiss_levels),4), size=(2000,1000))
    savefig(joinpath(outputs_dir,"$(loc_name).png"))
end