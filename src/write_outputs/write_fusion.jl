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
    dfFusion = inputs["dfFusion"]
    FUSION = inputs["FUSION"]

    T = inputs["T"]     # Number of time steps (hours)

    d = OrderedDict{String, Any}(
        "Imports" => vec(value.(EP[:vfusionimports][FUSION,1:T])),
        "Net Electric" => vec(value.(EP[:eFusionNetElec][FUSION,1:T])),
        "Gross Electric" => vec(value.(EP[:eTurbElec][FUSION,1:T])),
        "Recirc Power" => vec(value.(EP[:eRecircpwr][FUSION,1:T])),
        "Salt Heating" => vec(value.(EP[:vsaltpwr][FUSION,1:T])),
        "Fixed plant power" => vec(value.(EP[:eplantfix][FUSION,1:T])),
        "Var plant power" => vec(value.(EP[:eplantvar][FUSION,1:T])),
        "Reactor Thermal power" => vec(value.(EP[:vThermOutput][FUSION,1:T])),
        "Turbine Thermal input" => vec(value.(EP[:eTurbThermal][FUSION,1:T])),
        "Tritium inventory" => vec(value.(EP[:vtrit_inventory][FUSION,1:T])),
        "Tritium exports" => vec(value.(EP[:vtrit_exports][FUSION,1:T])),
        "Deuterium inventory" => vec(value.(EP[:vdeu_inventory][FUSION,1:T])),
        "Deuterium exports" => vec(value.(EP[:vdeu_exports][FUSION,1:T])),
        "Commitment State" => vec(value.(EP[:eFusionCommit][FUSION,1:T])),
        "Storage Inventory" => vec(value.(EP[:vThermStor][FUSION,1:T])),
        "Charge per Hour" => vec(value.(EP[:vThermChar][FUSION,1:T])),
        "Discharge per Hour" => vec(value.(EP[:vThermDis][FUSION,1:T])),
        "Net charge rate" => vec(value.(EP[:eThermStorNetDischarge][FUSION,1:T]))
    )

    fusion_df = DataFrame(d)    
    CSV.write(joinpath(path, "fusion_time.csv"), fusion_df)
end

function write_fusion_var(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    dfFusion = inputs["dfFusion"]
    FUSION = inputs["FUSION"]

    f_var = OrderedDict{String, Any}(
        "No. of Units" => value.(EP[:num_units][FUSION]),
        "Tritium Capacity" => value.(EP[:TritCap][FUSION]),
        "Deuterium Capacity" => value.(EP[:DeuCap][FUSION]),
        "Storage Capacity" => value.(EP[:vThermStorCap][FUSION]),
        "Discharge Capacity" => value.(EP[:vThermDisCap][FUSION])
    )

    fvar_df = DataFrame(f_var)
    CSV.write(joinpath(path, "fusion_var.csv"), fvar_df)
end

