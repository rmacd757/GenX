using CSV
using DataFrames
using XLSX
using OrderedCollections

current_dir = @__DIR__

dir_2_run = [
    # [joinpath(current_dir,"isone_1yr"), 0],
    # [joinpath(current_dir,"isone_1yr_simpleFusion"), 0],
    [joinpath(current_dir,"isone_1yr_noVRECap"), 0],
    [joinpath(current_dir,"isone_1yr_simpleFusion_noVRECap"), 0],
    [joinpath(current_dir,"isone_1yr_semiFusion_noVRECap"), 0],
]

df_2_save = OrderedDict{String,DataFrame}()

# for each dir in dir_2_run
# Check if dir/Results/cost.csv exists, if not, run the case
# if it does exist, skip the case
for x in dir_2_run
    dir = x[1]
    overwrite_flag = x[2]
    if !isfile(joinpath(dir,"Results","Primal_Capex_8500.0_EmissLevel_500.0","costs.csv")) || overwrite_flag == 1
        println("Running case in directory: ", dir)
        include(joinpath(dir,"Run.jl"))
    else
        println("Skipping case in directory: ", dir)
    end
end

for x in dir_2_run
    dir = x[1]
    # Load the costs.csv file
    cost_df = CSV.read(joinpath(dir,"Results","Primal_Capex_8500.0_EmissLevel_500.0","costs.csv"), DataFrame)
    # Load the capacity.csv file
    cap_df = CSV.read(joinpath(dir,"Results","Primal_Capex_8500.0_EmissLevel_500.0","capacity.csv"), DataFrame)
    # Load the capacityfactor.csv file
    capfac_df = CSV.read(joinpath(dir,"Results","Primal_Capex_8500.0_EmissLevel_500.0","capacityfactor.csv"), DataFrame)
    # Load the power.csv file
    power_df = CSV.read(joinpath(dir,"Results","Primal_Capex_8500.0_EmissLevel_500.0","power.csv"), DataFrame)
    
    temp = split(dir,"/")[end]
    case_name = join(split(temp,"_")[3:end],"_")

    df_2_save["$(case_name)_costs"] = cost_df
    df_2_save["$(case_name)_cap"] = cap_df
    df_2_save["$(case_name)_cap_fac"] = capfac_df
    df_2_save["$(case_name)_power"] = power_df
end

XLSX.writetable(joinpath(current_dir,"comparison-genx.xlsx"), overwrite=true, 
        df_2_save...
    )