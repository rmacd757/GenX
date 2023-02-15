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

    fusion_imports = vec(value.(EP[:vfusionimports][FUSION,1:T]))
    net_fusionpwr = vec(value.(EP[:eFusionPower][FUSION,1:T]))
    recirc_power = vec(value.(EP[:eRecircpwr][FUSION,1:T]))
    turbine_power = net_fusionpwr .+ recirc_power
    reactor_thermal = vec(value.(EP[:vThermOutput][FUSION,1:T]))

    tritium_inven = vec(value.(EP[:vtrit_inventory][FUSION,1:T]))
    deu_inven = vec(value.(EP[:vdeu_inventory][FUSION,1:T]))

    trit_exports = vec(value.(EP[:vtrit_exports][FUSION,1:T]))
    deu_exports = vec(value.(EP[:vdeu_exports][FUSION,1:T]))

    fixed_input = vec(value.(EP[:eplantfix][FUSION,1:T]))
    var_input = vec(value.(EP[:eplantvar][FUSION,1:T]))

    fusion_df = DataFrame(Time = 1:T, Imports = fusion_imports, Power = net_fusionpwr, Recirculate_Power = recirc_power, Turbine_Pwr = turbine_power, Reactor_Therm = reactor_thermal, Trit_Inven = tritium_inven, Deu_Inven = deu_inven, Trit_Exports = trit_exports, Deu_Exports = deu_exports, Fixed_Input = fixed_input, Var_Input = var_input)
    
    CSV.write(joinpath(path, "fusion.csv"), fusion_df)
end


