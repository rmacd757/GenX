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

    ## Initialize the variable for the reactor being on/standby (CAN probably find a way to call the inputs to get this value - most likely the COMMIT variable)
    @variable(EP, fusionstart[t=1:T], Bin)

    ## Create a variable for the fusion power nameplate capacity
    @variable(EP, NameplateCap)

    ## Convert nameplate capacity to effective thermal power capacity
    # Effective capacity is a constant fraction of the nameplate capacity
    # Formula is as follows: Effective Cap = NameplateCap*(Reactor Pulse Time/(Reactor Pulse Time + Dwell Time))
    # Dwell Time is 1 minute and reactor pulse time is 20 minutes
    @expression(EP, EffectiveCap, NameplateCap*(20/21))

    ## Reactor Thermal Output
    @variable(EP, ThermOutput[t=1:T])

    ## Constrain the Thermal Output with the Effective Capacity
    @constraint(EP, 0 <= ThermOutput[t]-EffectiveCap*fusionstart[t] <= 0)
end

### This function will calculate the recirculating power in the reactor
function recicrulate(EP::Model, inputs::Dict, setup::Dict)
    dfGen = inputs["dfGen"]

    T = inputs["T"]     # Number of time steps (hours)
    ## Define key variables
    @variable(EP, magcool == 10)  # Amount of power being delivered to cool the magnets 
    @expression(EP, plantfix[t=1:T], 10*fusionstart[t])  # Fixed power being delivered to the power plant (Number multiplied with the binary variable for start)
    @expression(EP, plantvar[t=1:T], 10*ThermOutput[t])  # Variable power being delivered to the power plant (Function of the ThermalOutput of the power plant)
    @variable(EP, saltpwr[t:1:T])  # Power being used to heat the salt 

    ## Combine to make recirculating power expression
    @expression(EP, Recircpwr[t=1:T], magcool+plantfix[t]+plantvar[t]+saltpwr[t])
end


### This function will calculate the thermal balance of the salt load_period_map
function saltthermal(EP::Model, inputs::Dict, setup::Dict)
    dfGen = inputs["dfGen"]

    T = inputs["T"]  # number of time steps (hours)

    ## Define key variables
    @variable(EP, salteff == 0.2) ## Salt electric heating efficiency
    @variable(EP, saltLosses[t=1:T])  ## Hourly thermal losses from salt loop  

    ## Expression for the thermal energy entering the turbine 
    @expression(EP, TurbThermal[t=1:T], saltpwr[t]*salteff+ThermOutput[t]-saltLosses[t])
end

