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


### This function will initialize the key variables in the fusion module 
function fusionthermalpower(EP::Model, inputs::Dict, setup::Dict)

    dfGen = inputs["dfGen"]

    T = inputs["T"]     # Number of time steps (hours)

    FUSION = inputs["FUSION"]
    COMMIT = inputs["COMMIT"] # For now, thermal resources are the only ones eligible for Unit Committment

	### Variables ###

    hr_unit = 3.412 ### hard coded mmbtu to mwh conversion

    con_eff = (dfGen[!,:Heat_Rate_MMBTU_per_MWh])./hr_unit   # Heat Rate of the technology is always greater than the mmbtu to mwh conversion
	## Decision variables for unit commitment
	# commitment state variable

    ## Expression to represent the commitment state of the fusion power generators
    @expression(EP, eFusionCommit[y in FUSION, t=1:T], EP[:vCOMMIT][y,t])

    # ## Initialize the variable for the reactor being on/standby (CAN probably find a way to call the inputs to get this value - most likely the COMMIT variable)
    # @variable(EP, fusioncommit[t=1:T], Bin)

    ## Fusion Power Nameplate Capacity (this is gross cap size)
    NameplateCap = (dfGen[!,:Cap_Size]).*con_eff   ## Convert gross capacity to thermal capacity (still keeping it in MW)

    ## Convert nameplate capacity to effective thermal power capacity
    # Effective capacity is a constant fraction of the nameplate capacity
    # Formula is as follows: Effective Cap = NameplateCap*(Reactor Pulse Time/(Reactor Pulse Time + Dwell Time))
    # Dwell Time is 1 minute and reactor pulse time is 20 minutes
    @expression(EP, eEffectiveCap[y in FUSION], NameplateCap[y]*(20/21))

    ## Reactor Thermal Output
    @variable(EP, vThermOutput[y in FUSION,t=1:T] >= 0)

    ## Constrain the Thermal Output with the Effective Capacity
    @constraint(EP, [y in FUSION,t=1:T], eEffectiveCap[y]*eFusionCommit[y,t] >= vThermOutput[y,t])
end

### This function will calculate the power provided by the reactor to the grid
function fusiongridpower(EP::Model, inputs::Dict, setup::Dict)
    T = inputs["T"]     # Number of time steps (hours)
    Z = inputs["Z"]     # Number of zones

    dfGen = inputs["dfGen"]
    dfFusion = inputs["dfFusion"]

    FUSION = inputs["FUSION"]
    
    hr_unit = 3.412 ### hard coded mmbtu to mwh conversion
    # con_eff = (dfGen[!,:Heat_Rate_MMBTU_per_MWh])./hr_unit   # Heat Rate of the technology is always greater than the mmbtu to mwh conversion

    ## Define key variables
    magcool = dfFusion[!,:Mag_Cool]  #Amount of power (MW) being delivered to cool the magnets

    @expression(EP, eplantfix[y in FUSION,t=1:T], 10*EP[:eFusionCommit][y,t])  # Fixed power being delivered to the power plant (Number multiplied with the binary variable for start)
    @expression(EP, eplantvar[y in FUSION,t=1:T], (1/12)*EP[:vThermOutput][y,t])  # Variable power being delivered to the power plant (Function of the ThermalOutput of the power plant)
    
    @variable(EP, vsaltpwr[y in FUSION,t=1:T] >= 0)  # Variable of the electric power being delivered to heat the salt

    ## Combine to make recirculating power expression
    @expression(EP, eRecircpwr[y in FUSION,t=1:T], magcool + eplantfix[y,t] + eplantvar[y,t] + vsaltpwr[y,t])

    ### Calculation for the thermal balance of the salt loop
    salteff = dfFusion[!,:Salt_Eff]      ## Salt electric heating efficiency (already converts MW to thermal)

    saltLosses = dfFusion[!,:Salt_Loss]     ## Hourly Thermal losses from salt loop (2 MW/hr)

    ## Expression for the thermal energy entering the turbine 
    @expression(EP, eTurbThermal[y in FUSION,t=1:T], vsaltpwr[y,t]*salteff + EP[:vThermOutput][y,t] - saltLosses)

    ## Define the turbine efficiency
    # turbeff = 0.40 #(Change out for heat rate in gen_data)

    ## Expression for the total amount of power delivered by the power plants to the grid
    @expression(EP, eFusionPower[y in FUSION,t=1:T], eTurbThermal[y,t]./((dfGen[y,:Heat_Rate_MMBTU_per_MWh])./hr_unit) - eRecircpwr[y,t])

    ## Variable for the amount of power that is being imported by the reactor when it is in standby mode
    @variable(EP, vfusionimports[y in FUSION,t=1:T] >= 0)

    ## Expression to sum the fusion imports across the zones 
    @expression(EP, efusionimports[t=1:T, z=1:Z], 
        sum(vfusionimports[y,t] for y in intersect(FUSION, dfGen[dfGen[!,:Zone].==z,:R_ID])))

    EP[:ePowerBalance] -= efusionimports

    ## Constraints for the electric power node
    @constraint(EP, [y in FUSION, t=1:T], eFusionPower[y,t] == EP[:vP][y,t])
    @constraint(EP, PwrNode[y in FUSION, t=1:T], eFusionPower[y,t] + vfusionimports[y,t] == EP[:vP][y,t] + eRecircpwr[y,t])
end


function fusionfuel(EP::Model, inputs::Dict, setup::Dict)
    dfGen = inputs["dfGen"]

    T = inputs["T"]     # Number of time steps (hours)

    FUSION = inputs["FUSION"]
    dfFusion = inputs["dfFusion"]

    ##### Tritium Balance #####

    trit_stor_cap = dfFusion[!,:Tritium_Cap]

    # Tritium inventory
    @variable(EP, trit_inventory[y in FUSION, t=0:T], lower_bound == 0, upper_bound == trit_stor_cap)
  
    # Fuel Coefficient, Breeding Coefficient, Loss Rate
    trit_fuel = dfFusion[!,:Trit_Fuel]
    trit_breed = dfFusion[!,:Trit_Breed]
    trit_loss = dfFusion[!,:Trit_Loss]

    # Tritium Balance
    @expression(EP, trit_balance[y in FUSION, t=1:T], trit_inventory[y,t-1] + EP[:vThermOutput][y,t]*(trit_breed[y,t] - trit_fuel[y,t]) - trit_loss[y,t])
    

    ##### Deuterium Balance #####

    deu_stor_cap = dfFusion[!,:Deu_Cap]

    # Deuterium inventory
    @variable(EP, deu_inventory[y in FUSION, t=0:T], lower_bound == 0, upper_bound == deu_stor_cap)

    # Fuel Coefficient, Loss Rate
    deu_fuel = dfFusion[!,:Deu_Fuel]
    deu_loss = dfFusion[!,:Deu_Loss]

    # Deuterium exports
    @variable(EP, deu_exports[y in FUSION, t=1:T])

    # Deuterium balance
    @expression(EP, deu_balance[y in FUSION, t=1:T], due_inventory[y,t-1] - EP[:vThermOutput][y,t]*(deu_fuel[y,t]) - deu_loss[y,t] - deu_exports[y,t])
end

## This function is the overall function for fusion power 
function fusion!(EP::Model, inputs::Dict, setup::Dict)
    fusionthermalpower(EP,inputs,setup)
    fusiongridpower(EP,inputs,setup)
    fusionfuel(EP,inputs,setup)
end