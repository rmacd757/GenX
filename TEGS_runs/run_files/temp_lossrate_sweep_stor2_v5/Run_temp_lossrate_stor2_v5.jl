using GenX
using CSV
using DataFrames
using JuMP

function tegs_case(dsratio, outputs_path_case, mysetup, myinputs, OPTIMIZER)
    # Calculate and save baseline emissions
    mkpath(outputs_path_case)
    outputs_path = joinpath(outputs_path_case, "$(dsratio)")
    
    println("Generating the Optimization Model")
    EP = generate_model(mysetup, myinputs, OPTIMIZER)

    # {Discharge capacity} * {dsratio} == {energy capacity}
    @constraint(EP, TEGSdsratio, EP[:eTotalCap][end] * dsratio == EP[:eTotalCapEnergy][end])

    # {Discharge power} / {discharge capacity} <= {energy} / {energy capacity}
    # {Discharge power} / {discharge capacity} <= {energy} / {discharge capacity * dsratio}
    # {Discharge power} * {dsratio} <= {energy}
    @constraint(EP, TEGSeffpower[t=1:myinputs["T"]], EP[:vP][end,t] * dsratio <= EP[:vS][end,t])

    println("Solving Model")
    EP, solve_time = solve_model(EP, mysetup)
    myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

    # Run MGA if the MGA flag is set to 1 else only save the least cost solution
    println("Writing Output")
    elapsed_time = @elapsed write_outputs(EP, outputs_path, mysetup, myinputs)
    println("Time elapsed for writing is")
    println(elapsed_time)

    return objective_value(EP)

end

function tegs_brent(x0::Number, x1::Number, outputs_path_case::String, mysetup::Dict{Any, Any}, myinputs::Dict{Any, Any}, OPTIMIZER::MOI.OptimizerWithAttributes, xtol::AbstractFloat=1e-7, ytol=2eps(Float64), maxiter::Integer=50)
    EPS = eps(Float64)
    y0 = tegs_case(x0, outputs_path_case, mysetup, myinputs, OPTIMIZER)
    y1 = tegs_case(x1, outputs_path_case, mysetup, myinputs, OPTIMIZER)
    if abs(y0) < abs(y1)
        # Swap lower and upper bounds.
        x0, x1 = x1, x0
        y0, y1 = y1, y0
    end
    x2 = x0
    y2 = y0
    x3 = x2
    bisection = true
    for _ in 1:maxiter
        # x-tolerance.
        if abs(x1-x0) < xtol
            return x1
        end

        # Use inverse quadratic interpolation if f(x0)!=f(x1)!=f(x2)
        # and linear interpolation (secant method) otherwise.
        if abs(y0-y2) > ytol && abs(y1-y2) > ytol
            x = x0*y1*y2/((y0-y1)*(y0-y2)) +
            x1*y0*y2/((y1-y0)*(y1-y2)) +
            x2*y0*y1/((y2-y0)*(y2-y1))
        else
            x = x1 - y1 * (x1-x0)/(y1-y0)
        end

        # Use bisection method if satisfies the conditions.
        delta = abs(2EPS*abs(x1))
        min1 = abs(x-x1)
        min2 = abs(x1-x2)
        min3 = abs(x2-x3)
        if (x < (3x0+x1)/4 && x > x1) ||
            (bisection && min1 >= min2/2) ||
            (!bisection && min1 >= min3/2) ||
            (bisection && min2 < delta) ||
            (!bisection && min3 < delta)
            x = (x0+x1)/2
            bisection = true
        else
            bisection = false
        end

        y = tegs_case(x, outputs_path_case, mysetup, myinputs, OPTIMIZER)
        # y-tolerance.
        if abs(y) < ytol
            return x
        end
        x3 = x2
        x2 = x1
        if sign(y0) != sign(y)
            x1 = x
            y1 = y
        else
            x0 = x
            y0 = y
        end
        if abs(y0) < abs(y1)
            # Swap lower and upper bounds.
            x0, x1 = x1, x0
            y0, y1 = y1, y0
        end
    end
    error("Max iteration exceeded")
end

function ITP_root(f, a, b, ϵ=eps(), κ₁=0.2/(b-a), κ₂=2, n₀=1)
    0 < κ₁ < Inf             ||   error("κ₁ must be between 0 and ∞")
    1 ≤ κ₂ < 1 + (1 + √5)/2   ||   error("κ₂ must be between 1 and 1+ (1 + √5)/2 (1 + the golden ratio)")
    0 < n₀ < Inf             ||   error("n₀ must be between 0 and ∞")
    n_1div2 = ceil(Int, log2((b-a)/2ϵ))
    nₘₐₓ = n_1div2 + n₀
    y_a = f(a)
    y_b = f(b)
    sign(y_a) == sign(y_b)  &&  error("sign(f(a)) = sign(f(b)). There is no guaranteed root in the given interval.")
    j = 0
    while b-a > 2ϵ
        # Calculating parameters:
        x_1div2 = (a+b)/2
        r = ϵ*2^(nₘₐₓ - j) - (b-a)/2
        δ = κ₁*(b-a)^κ₂
        
        # Interpolation:
        x_f = /(y_b*a - y_a*b, y_b-y_a)
        
        # Truncation:
        σ = sign(x_1div2 - x_f)
        δ ≤ abs(x_1div2 - x_f) ? (x_t=x_f+σ*δ) : (x_t = x_1div2)
        
        # Projection:
        abs(x_t - x_1div2) ≤ r ? (x_ITP = x_t) : (x_ITP = x_1div2 - σ*r)
        
        # Updating Interval:
        y_ITP = f(x_ITP)
        if y_ITP > 0
            b = x_ITP
            y_b = y_ITP
        elseif y_ITP < 0
            a = x_ITP
            y_a = y_ITP
        else
            a = b = x_ITP
        end
        j += 1
    end
    return (a+b)/2
end

function ITP_min(f, a, b, args, ϵ=eps(), κ₁=0.2/(b-a), κ₂=2, n₀=1)
    0 < κ₁ < Inf             ||   error("κ₁ must be between 0 and ∞")
    1 ≤ κ₂ < 1 + (1 + √5)/2   ||   error("κ₂ must be between 1 and 1+ (1 + √5)/2 (1 + the golden ratio)")
    0 < n₀ < Inf             ||   error("n₀ must be between 0 and ∞")

    n_1div2 = ceil(Int, log2((b-a)/2ϵ))
    nₘₐₓ = n_1div2 + n₀

    c = (a+b) / 2
    y_a = f(a, args...)
    y_c = f(c, args...)
    y_b = f(b, args...)
    # sign(y_a) == sign(y_b)  &&  error("sign(f(a)) = sign(f(b)). There is no guaranteed root in the given interval.")

    println("A = $(a)")
    println("C = $(b)")
    println("B = $(c)")

    j = 0
    while b-a > 2ϵ

        # Calculating parameters:
        x_1div2 = (a+b)/2
        r = ϵ*2^(nₘₐₓ - j) - (b-a)/2
        δ = κ₁*(b-a)^κ₂
        
        # Interpolation:
        x_f = /(y_b*a - y_a*b, y_b-y_a)
        
        # Truncation:
        σ = sign(x_1div2 - x_f)
        δ ≤ abs(x_1div2 - x_f) ? (x_t=x_f - σ * δ) : (x_t = x_1div2)
        
        # Projection:
        abs(x_t - x_1div2) ≤ r ? (x_ITP = x_t) : (x_ITP = x_1div2 + σ * r)

        println("A = $(a)")
        println("C = $(b)")
        println("B = $(c)")
        println("X = $(x_ITP)")
        
        # Updating Interval:
        y_ITP = f(x_ITP, args...)
        if y_c < y_ITP
            # (a,c,b) -> (a,c,x)
            b = x_ITP
            y_b = y_ITP
        elseif y_ITP < y_c
            # (a,c,b) -> (c,x,b)
            a = c
            y_a = y_c
            c = x_ITP
            y_c = y_ITP
        else
            a = b = x_ITP
        end
        j += 1

    end
    # return (a+b)/2
    return c
end

############################################
# Some necessary helper functions
############################################
include("../RunTools.jl")

function costpermwdischarge(T)
    T = Int(T)
    discharge_vs_temp = Dict{Float64, Float64}(
        2400 => 0.3097 * 1e6,
        2300 => 0.3427 * 1e6,
        2100 => 0.4926 * 1e6,
        1900 => 0.8924 * 1e6,
    )
    return discharge_vs_temp[T]
end

function costpermwcharge(T)
    return 0.0185 * 1e6
end

function costpermwhstor(T, lossrate)
    storage_vs_temp = Dict{Float64, Array{Float64}}(
        2400 => 1e3 * [15.827, 12.0636, 10.8529, 10.2619, 9.9158, 9.6914, 9.5363, 9.4243, 9.3411, 9.2779, 9.2292, 9.1915, 9.1622, 9.1395, 9.122],
        2300 => 1e3 * [15.1493, 11.7375, 10.6382, 10.1019, 9.7883, 9.5854, 9.4456, 9.3451, 9.2707, 9.2146, 9.1717, 9.1388, 9.1136, 9.0943, 9.0799],
        2100 => 1e3 * [13.8087, 11.0887, 10.2103, 9.7827, 9.5337, 9.3737, 9.2644, 9.1867, 9.1301, 9.0881, 9.0568, 9.0336, 9.0165, 9.0042, 8.9958],
        1900 => 1e3 * [12.4875, 10.4445, 9.7845, 9.4646, 9.2799, 9.1626, 9.0836, 9.0287, 8.9897, 8.9618, 8.9421, 8.9284, 8.9194, 8.9141, 8.9117],
    )
    return storage_vs_temp[T][floor(Int, lossrate)]
end

function annuitize(cost::Float64, lifetime::Float64, discount_rate::Float64)
    return cost * discount_rate / (1 - (1 + discount_rate)^(-lifetime))
end

function TEGS_costs(stortype::Int, T::Float64, lossrate::Float64, lifetime::Float64, discount_rate::Float64)
    discharge_cost = annuitize(costpermwdischarge(T), lifetime, discount_rate)
    storage_cost = annuitize(costpermwhstor(T, lossrate), lifetime, discount_rate)
    charge_cost = annuitize(costpermwcharge(T), lifetime, discount_rate)
    if stortype == 1
        costs = Dict{String, Float64}(
        "Inv_Cost_per_MWyr" => discharge_cost + charge_cost,
        "Fixed_OM_Cost_per_MWyr" => 0.25 * (discharge_cost + charge_cost),
        "Inv_Cost_per_MWhyr" => storage_cost,
        "Fixed_OM_Cost_per_MWhyr" => 0.25 * storage_cost,
    )
    elseif stortype == 2
        costs = Dict{String, Float64}(
        "Inv_Cost_per_MWyr" => discharge_cost,
        "Fixed_OM_Cost_per_MWyr" => 0.25 * discharge_cost,
        "Inv_Cost_per_MWhyr" => storage_cost,
        "Fixed_OM_Cost_per_MWhyr" => 0.25 * storage_cost,
        "Inv_Cost_Charge_per_MWyr" => charge_cost,
        "Fixed_OM_Cost_Charge_per_MWyr" => 0.25 * charge_cost,
    )
    else
        error("Invalid storage type")
    end
    return costs
end

function selectresource(dfGen::DataFrame, resourcename::String)
    dfgen_search = groupby(dfGen, :Resource)
    return dfgen_search[(resourcename,)]
end
    
function setTEGScosts!(dfGen::DataFrame, stortype::Int, T::Float64, lossrate::Float64, lifetime::Float64, discount_rate::Float64)
    TEGS_input = selectresource(dfGen, "TEGS")
    tegs_costs = TEGS_costs(stortype, T, lossrate, lifetime, discount_rate)
    for (param, value) in tegs_costs
        TEGS_input[!,param] .= value
    end   
end

############################################
# Case Definitions 
# All cases intended to be run from the run-file directory
############################################
root_dir = dirname(dirname(dirname(@__FILE__))) # Should be ../TEGS_runs
run_name = "temp_lossrate_sweep_stor2_v5"
dropbox_path = "/Users/rmacd/Dropbox/1_Academics/Research/22-TEGS_modelling/TEGS GenX shared folder/GenX_runs"
# dropbox_path = "D:/Dropbox/1_Academics/Research/22-TEGS_modelling/TEGS GenX shared folder/GenX_runs"
# dropbox_path = "/media/rmacd/LargeHD/Dropbox/1_Academics/Research/22-TEGS_modelling/TEGS GenX shared folder/GenX_runs"

location_dir = Dict{String, String}(
    "newEngland" => joinpath(root_dir, "data", "newEngland_stor2_v2"),
    "texas" => joinpath(root_dir, "data", "texas_stor2"),
)

logging_notes = Array{String, 1}()

emissions_levels = Dict{String, Float64}(
    "0.20"=>0.2,
    "0.10"=>0.1, 
    "0.05"=>0.05,
    "0.01"=>0.01
    )

emiss_targets = Dict{String, Dict{String, Float64}}(
    "newEngland" => Dict{String, Float64}(),
    "texas" => Dict{String, Float64}(),
)

for (loc_name, loc_emiss) in emiss_targets
    outputs_path = joinpath(dropbox_path, "outputs", loc_name, "emissions_and_baseline_v2", "baseline")
    # outputs_path = joinpath(root_dir, "outputs", loc_name, "emissions_and_baseline", "baseline")
    emissions_df = DataFrame(CSV.File(joinpath(outputs_path, "emissions.csv")))
    baseline_emiss = emissions_df.Total[1]
    for (emiss_name, emiss_level) in emissions_levels
        loc_emiss[emiss_name] = baseline_emiss * emiss_level
    end
end
        
temperatures = Array{Float64}([2400, 2300, 2100, 1900])
lossrates = Array{Float64}([3, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])

for (loc_name, loc_path) in location_dir
    case = loc_path
    genx_settings = get_settings_path(case, "genx_settings.yml") #Settings YAML file path
    mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters

    mysetup["CO2Cap"] = 1
    push!(logging_notes, "Emissions cap activated\n")

    scale_factor = mysetup["ParameterScale"] == 1 ? ModelScalingFactor : 1

    inputs_path = case
    settings_path = get_settings_path(case)

    # outputs_path_root = joinpath(root_dir, "outputs", loc_name, run_name)
    outputs_path_root = joinpath(dropbox_path, "outputs", loc_name, run_name)
    if !isdir(outputs_path_root)
        mkpath(outputs_path_root)
    end

    ### Cluster time series inputs if necessary and if specified by the user
    TDRpath = joinpath(case, mysetup["TimeDomainReductionFolder"])

    if mysetup["TimeDomainReduction"] == 1
        if !time_domain_reduced_files_exist(TDRpath)
            println("Clustering Time Series Data (Grouped)...")
            cluster_inputs(inputs_path, settings_path, mysetup)
        else
            println("Time Series Data Already Clustered.")
        end
    end

    ### Configure solver
    println("Configuring Solver")
    OPTIMIZER = configure_solver(mysetup["Solver"], settings_path)

    ### Load inputs
    println("Loading Inputs")
    myinputs = load_inputs(mysetup, case)

    LIFETIME = 30.
    DISCOUNT_RATE = 0.05
    STOR_TYPE = 2

    push!(logging_notes, "LIFETIME = $LIFETIME\n")
    push!(logging_notes, "DISCOUNT_RATE = $DISCOUNT_RATE\n")
    push!(logging_notes, "STOR_TYPE = $STOR_TYPE\n")
    
    #### Running a case
    for (emiss_name, emiss_target) in emiss_targets[loc_name]
        myinputs["dfMaxCO2"] = emiss_target * 1e3 / scale_factor
        for T in temperatures
            for lossrate in lossrates
                push!(logging_notes, "Running $(emiss_name)_$(T)_$(lossrate)_stor2 case\n")
                setTEGScosts!(myinputs["dfGen"], STOR_TYPE, T, lossrate, LIFETIME, DISCOUNT_RATE)
                TEGS_input = selectresource(myinputs["dfGen"], "TEGS")
                TEGS_input[!, "Self_Disch"] .= lossrate / 100 / 24. # Convert daily loss rate to hourly
                TEGS_input[!, "STOR"] .= 2

                # Calculate and save baseline emissions
                outputs_path = joinpath(outputs_path_root, "$(emiss_name)_$(T)_$(lossrate)_stor2")
                if isdir(outputs_path)
                    println("$(emiss_name)_$(T)_$(lossrate)_stor2 already exists. Skipping...")
                    continue
                end

                # opt_dsratio = tegs_brent(1.0, 150.0, outputs_path, mysetup, myinputs, OPTIMIZER, 10.0)
                opt_dsratio = ITP_min(tegs_case, 1.0, 100.0, (outputs_path, mysetup, myinputs, OPTIMIZER), 10.0)

                push!(logging_notes, "1 : $(opt_dsratio), discharge : storage capacity was optimal. Saving as $(string(outputs_path, "_ds$(opt_dsratio)"))\n")

                cp(joinpath(outputs_path, "$(opt_dsratio)"), string(outputs_path, "_ds$(dsratio)"))

            end
        end
    end
end

# Write logging notes
open(joinpath(dirname(@__FILE__), "$(run_name).log"), "w") do f
    for line in logging_notes
        write(f, line)
    end
end