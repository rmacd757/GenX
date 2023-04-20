"""
GenX: An Configurable Capacity Expansion Model
Copyright (C) 2021,  Massachusetts Institute of Technology
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
A complete copy of the GNU General Public License v2 (GPLv2) is available
in LICENSE.txt.  Users uncompressing this from an archive may not have
received this license file.  If not, see <http://www.gnu.org/licenses/>.
"""

function write_fusion(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    # dfFusion = inputs["dfFusion"]
    FUSION = inputs["FUSION"]
    T = inputs["T"]     # Number of time steps (hours)

    # Make ~/fusion directory if it doesn't exist
    fusion_results_path = joinpath(path, "fusion")
    if !isdir(fusion_results_path)
        mkpath(fusion_results_path)
    end

    for y in FUSION
        result_filename = string("fusion_time_", y, ".csv")
        fusion_timeseries = OrderedDict{String, Any}(
            "Imports" => :vfusionimports,
            "Net Electric" => :eFusionNetElec,
            "Gross Electric" => :eTurbElec,
            "Recirc Power" => :eRecircpwr,
            "Salt Heating" => :vsaltpwr,
            "Fixed plant power" => :eplantfix,
            "Var plant power" => :eplantvar,
            "Reactor Thermal power" => :vThermOutput,
            "Turbine Thermal input" => :eTurbThermal,
            "Tritium inventory" => :vtrit_inventory,
            "Tritium exports" => :vtrit_exports,
            "Deuterium inventory" => :vdeu_inventory,
            "Deuterium imports" => :vdeu_imports,
            "Commitment State" => :eFusionCommit,
            "Storage Inventory" => :vThermStor,
            "Charge per Hour" => :vThermChar,
            "Discharge per Hour" => :vThermDis,
            "Net discharge rate" => :eThermStorNetDischarge
        )
        for (name, model_key) in fusion_timeseries
            if haskey(EP, model_key)
                fusion_timeseries[name] = get_fusion_data_vector(EP[model_key], y, T)
            else
                fusion_timeseries[name] = zeros(T)
            end
        end
        fusion_df = DataFrame(fusion_timeseries)    
        CSV.write(joinpath(fusion_results_path, result_filename), fusion_df)
    end
end

function get_fusion_data_vector(model_data::JuMP.Containers.DenseAxisArray{AffExpr, 2, Tuple{Vector{Int64}, Base.OneTo{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}, JuMP.Containers._AxisLookup{Base.OneTo{Int64}}}}, R_ID::Int64, maxtime::Int64, mintime::Int64=1)
    return value.(model_data[R_ID, mintime:maxtime]).data
end

function get_fusion_data_vector(model_data::JuMP.Containers.DenseAxisArray{VariableRef, 2, Tuple{Vector{Int64}, Base.OneTo{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}, JuMP.Containers._AxisLookup{Base.OneTo{Int64}}}}, R_ID::Int64, maxtime::Int64, mintime::Int64=1)
    return value.(model_data[R_ID, mintime:maxtime]).data
end

function write_fusion_var(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    dfFusion = inputs["dfFusion"]
    FUSION = inputs["FUSION"]
    dfGen = inputs["dfGen"]

    numunitf = zeros(size(inputs["dfGen"]))
    reactorcap = zeros(size(inputs["dfGen"]))
    turbcap = zeros(size(inputs["dfGen"]))
    tritcap = zeros(size(inputs["dfGen"]))
    deucap = zeros(size(inputs["dfGen"]))
    fthermstorcap = zeros(size(inputs["dfGen"]))
    fthermdiscap = zeros(size(inputs["dfGen"]))
    
    for i in FUSION
        numunitf[i] = value(EP[:num_units][i])
        reactorcap[i] = value(EP[:eTotalCap][i])
        turbcap[i] = value(EP[:vTurbElecCap][i])
        tritcap[i] = value(EP[:vTritCap][i])
        deucap[i] = value(EP[:vDeuCap][i])  
        if any(dfFusion[!,:Add_Therm_Stor].>0)
            fthermstorcap[i] = value(EP[:vThermStorCap][i])
            fthermdiscap[i] = value(EP[:vThermDisCap][i])
        else
            fthermstorcap[i] = 0
            fthermdiscap[i] = 0
        end
    end

    f_var = DataFrame(
        Fusion_tech = inputs["FUSION"],
        Num_Units = numunitf[FUSION],
        Reactor_Gross_Cap = reactorcap[FUSION],
        Turbine_Cap = turbcap[FUSION],
        Trit_Cap = tritcap[FUSION],
        Deu_Cap = deucap[FUSION],
        Therm_Stor_Cap = fthermstorcap[FUSION],
        Therm_Dis_Cap = fthermdiscap[FUSION],
    )

    CSV.write(joinpath(path, "fusion_var.csv"), f_var)
end
