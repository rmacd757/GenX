using GenX
using JuMP

using OrderedCollections
using DataFrames
using CSV

# case = dirname(@__FILE__)
case = joinpath(pwd(),"fusion_runs","data","primal_6zoneandQC_baseline")

function get_settings_path(case::AbstractString)
    return joinpath(case, "Settings")
end

function get_settings_path(case::AbstractString, filename::AbstractString)
    return joinpath(get_settings_path(case), filename)
end

function get_default_output_folder(case::AbstractString)
    return joinpath(case, "Results")
end

genx_settings = get_settings_path(case, "genx_settings.yml") #Settings YAML file path
mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters

inputs_path = case
settings_path = get_settings_path(case)

### Cluster time series inputs if necessary and if specified by the user
TDRpath = joinpath(case, mysetup["TimeDomainReductionFolder"])

function time_domain_reduced_files_exist(tdrpath)
    tdr_load = isfile(joinpath(tdrpath,"Load_data.csv"))
    tdr_genvar = isfile(joinpath(tdrpath,"Generators_variability.csv"))
    tdr_fuels = isfile(joinpath(tdrpath,"Fuels_data.csv"))
    return (tdr_load && tdr_genvar && tdr_fuels)
end

if mysetup["TimeDomainReduction"] == 1
    if !time_domain_reduced_files_exist(TDRpath)
        println("Clustering Time Series Data (Grouped)...")
        cluster_inputs(inputs_path, settings_path, mysetup)
    else
        println("Time Series Data Already Clustered.")
    end
end

function write_fusion_duals(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    f_var = DataFrame(
        Gross_Cap_Dual = -dual(EP[:fusiongrosseleccap]),
        # Turbine_Dual = dual(EP[:fusionturbeleccap]),
        # Stor_Dual = dual(EP[:fusionthermalstorcap]),
        # Stor_Dis_Dual = dual(EP[:fusionthermaldiscap])
    )
    CSV.write(joinpath(path, "fusion_duals.csv"), f_var)
end

### Configure solver
println("Configuring Solver")
OPTIMIZER = configure_solver(mysetup["Solver"], settings_path)
set_optimizer_attribute(OPTIMIZER, "BarHomogeneous", 1)

#### Running a case

### Load inputs
println("Loading Inputs")
myinputs = load_inputs(mysetup, case)

# reactor_cap_list = 1000. .* [20, 17.5, 15, 12.5, 10., 7.5, 5., 2.5, 1.]
# turb_ratio_list = 1.0 .* [4, 3, 2, 1]
reactor_cap_list = range(1000., stop=20000., step=1000.)
turb_ratio_list = [1.0]
# turb_ratio_list = range(1., stop=4., step=0.5)
therm_stor_list = [1500.0]
therm_dis_ratio_list = [300.0]
emiss_lim_list = 100.0 .* [2.5, 5, 7.5, 10, 15, 20, 25]

# reactor_cap = 1000. * 20.
# turb_ratio = 1.0 * 4.
# therm_stor = 1500.0
# therm_dis_ratio = 300.0

mysetup["CO2Cap"] = 1
scale_factor = mysetup["ParameterScale"] == 1 ? ModelScalingFactor : 1

myinputs["dfMaxCO2"] = emiss_lim * 1e3 / scale_factor
outputs_path = joinpath(case, "Results", "EmissLevel_" * string(emiss_lim) * "_ReactorCap_" * string(reactor_cap) * "_TurbRatio_" * string(turb_ratio) * "_ThermStor_" * string(therm_stor) * "_ThermDisRatio_" * string(therm_dis_ratio))
# elapsed_time = @elapsed write_outputs(EP, outputs_path, mysetup, myinputs)

mkpath(outputs_path)

println("Generating the Optimization Model")
EP = generate_model(mysetup, myinputs, OPTIMIZER)

println("Solving Model")
EP, solve_time = solve_model(EP, mysetup)
myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

# Run MGA if the MGA flag is set to 1 else only save the least cost solution
println("Writing Output")
outputs_path = get_default_output_folder(case)