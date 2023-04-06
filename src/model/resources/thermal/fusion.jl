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
    dfFusion = inputs["dfFusion"]

    T = inputs["T"]     # Number of time steps (hours)

    FUSION = inputs["FUSION"]

	### Variables ###
    MMBTU_PER_MWH = 3.412 ### hard coded mmbtu to mwh conversion

    con_mWh_heatrt = dfGen[!,:Heat_Rate_MMBTU_per_MWh] ./ MMBTU_PER_MWH   # Heat Rate of the technology is always greater than the mmbtu to mwh conversion

    # Expression to represent the commitment state of the fusion power generators
    # Note, this is the number of commited units in the cluster if using linearized unit commitment
    @expression(EP, eFusionCommit[y in FUSION, t=1:T], EP[:vCOMMIT][y,t])

    ## Convert nameplate capacity to effective thermal power capacity
    # Effective capacity is a constant fraction of the nameplate capacity
    # Formula is as follows: Effective Cap = NameplateCap*(Reactor Pulse Time/(Reactor Pulse Time + Dwell Time))
    # Dwell Time is 1 minute and reactor pulse time is 20 minutes

    ## Reactor Thermal Output
    @variable(EP, vThermOutput[y in FUSION,t=1:T] >= 0)

    ## Constrain the Thermal Output with the Effective Capacity
    pulse_ratio = fusionpulseratio.(dfFusion[!,:Pulse_length_min], dfFusion[!,:Pulse_gap_min])
    @constraint(EP, [y in FUSION,t=1:T], vThermOutput[y,t] <= eFusionCommit[y,t] * dfGen[!,:Cap_Size][y] * con_mWh_heatrt[y] * pulse_ratio[y])
end

### This function will calculate the power provided by the reactor to the grid
function fusiongridpower(EP::Model, inputs::Dict, setup::Dict)
    T = inputs["T"]     # Number of time steps (hours)
    Z = inputs["Z"]     # Number of zones

    dfGen = inputs["dfGen"]
    dfFusion = inputs["dfFusion"]

    FUSION = inputs["FUSION"]
    
    # This calculation finds the number of units based on the gross power of the reactor, not net
    @expression(EP, num_units[y in FUSION], EP[:eTotalCap][y] / dfGen[!,:Cap_Size][y])

    MMBTU_PER_MWH = 3.412 ### hard coded mmbtu to mwh conversion

    ## Turbine efficiency
    turbefficiency = MMBTU_PER_MWH ./ dfGen[!,:Heat_Rate_MMBTU_per_MWh]

    ## Define key variables
    magcool = dfFusion[!,:Mag_Cool]  #Amount of power (MW) being delivered to cool the magnets

    ### Calculation for the thermal balance of the salt loop
    salteff = dfFusion[!,:Salt_Eff]      ## Salt electric heating efficiency (already converts MW to thermal)
    saltLosses = dfFusion[!,:Salt_Loss]  ## Hourly Thermal losses from salt loop (2 MW/hr) for each reactor

    # 10 MW is a good guess for fix plant power
    @expression(EP, eplantfix[y in FUSION,t=1:T], dfFusion[y,:Plant_fixed_power_MWe] * EP[:eFusionCommit][y,t])  # Fixed power being delivered to the power plant (Number multiplied with the binary variable for start)
    @expression(EP, eplantvar[y in FUSION,t=1:T], dfFusion[y,:Plant_var_power_MWeperMWt] * EP[:vThermOutput][y,t])  # Variable power being delivered to the power plant (Function of the ThermalOutput of the power plant)
    
    @variable(EP, 0 <= vsaltpwr[y in FUSION,t=1:T])  # Variable of the electric power being delivered to heat the salt
    @constraint(EP, [y in FUSION,t=1:T], vsaltpwr[y,t] <= num_units[y] * saltLosses[y] / salteff[y])

    ## Combine to make recirculating power expression
    @expression(EP, eRecircpwr[y in FUSION,t=1:T], num_units[y] * magcool[y] + eplantfix[y,t] + vsaltpwr[y,t] + eplantvar[y,t])
    @constraint(EP, eRecircpwr .>= 0)

    ## Expression for the thermal energy entering the turbine 
    # @expression(EP, eTurbThermal[y in FUSION,t=1:T], vsaltpwr[y,t] * salteff[y] + EP[:vThermOutput][y,t] - num_units[y] * saltLosses[y])
    # @constraint(EP, eTurbThermal .>= 0)
    # @variable(EP, eTurbThermal[y in FUSION,t=1:T] >= 0.0)
    # @constraint(EP, [y in FUSION,t=1:T], eTurbThermal[y,t] == EP[:eThermStorNetDischarge][y,t] + EP[:vThermOutput][y,t] - num_units[y] * saltLosses[y] + EP[:vsaltpwr][y,t] * salteff[y])
    @expression(EP, eTurbThermal[y in FUSION,t=1:T], EP[:eThermStorNetDischarge][y,t] + EP[:vThermOutput][y,t] - num_units[y] * saltLosses[y] + EP[:vsaltpwr][y,t] * salteff[y])
    @constraint(EP, eTurbThermal .>= 0)

    ## Expression for the amount of electric power coming out from the turbine
    # First, define the variable turbine gross electric power capacity
    # @variable(EP, vTurbElecCap[y in FUSION], lower_bound = 0.)

    # Then, define the expression for the amount of electric power exiting turbine
    @expression(EP, eTurbElec[y in FUSION,t=1:T], eTurbThermal[y,t] .* turbefficiency[y])


    # Place constraints on the amount of energy exiting the turbine
    # @constraint(EP, [y in FUSION,t=1:T], eTurbElec[y,t] <= vTurbElecCap[y])
    # @constraint(EP, eTurbElec[y,t] .>= 0)

    # # Place constraint on turbine minimum power
    # @constraint(EP, [y in FUSION,t=1:T], eTurbElec[y,t] <= vTurbElecCap[y,t] * EP[:eFusionCommit][y,t])
    # @constraint(EP, [y in FUSION,t=1:T], eTurbElec[y,t] >= 0.4 * vTurbElecCap[y,t] * EP[:eFusionCommit][y,t])

    # # Then, add the costs of the turbine to the fixed costs
    # turb_cost = dfFusion[!,:Turb_Cap_Cost]

    # # Adding costs of variable turbine capacity to fixed costs
    # EP[:eCFix][FUSION] .+= (vTurbElecCap[FUSION].*turb_cost[FUSION])
    # @expression(EP, eSumTurbFix, (sum((vTurbElecCap[y].*turb_cost[y] for y in FUSION))))
    # EP[:eTotalCFix] += eSumTurbFix
    # EP[:eObj] += eSumTurbFix 

    # Adding costs of fixed turbine capacity to fixed costs
    # EP[:eCFix][FUSION] .+= (eMaxElec[FUSION].*turb_cost[FUSION])
    # @expression(EP, eSumTurbFix, (sum((eTurbElec[y].*turb_cost[y] for y in FUSION))))
    # EP[:eTotalCFix] += eSumTurbFix
    # EP[:eObj] += eSumTurbFix 

    ## Expression for the amount of net power that is actually going to the grid
    @expression(EP, eFusionNetElec[y in FUSION,t=1:T], eTurbElec[y,t] - eRecircpwr[y,t])
    @constraint(EP, eFusionNetElec .>= 0)

    ## Variable for the amount of power that is being imported by the reactor when it is in standby mode
    @variable(EP, vfusionimports[y in FUSION,t=1:T] >= 0)
    @constraint(EP, [y in FUSION,t=1:T], vfusionimports[y,t] <= eRecircpwr[y,t])

    ## Expression to sum the fusion imports across the zones 
    @expression(EP, efusionimports[t=1:T, z=1:Z], 
        sum(vfusionimports[y,t] for y in intersect(FUSION, dfGen[dfGen[!,:Zone].==z,:R_ID])))

    EP[:ePowerBalance] -= efusionimports

    ## Constraints for the electric power node (Supply of power on the left, Demand for power on right)
    @constraint(EP, PwrNode[y in FUSION, t=1:T], eTurbElec[y,t] + vfusionimports[y,t] == EP[:vP][y,t] + EP[:eRecircpwr][y,t])
end


function fusionfuel(EP::Model, inputs::Dict, setup::Dict)
    T = inputs["T"]     # Number of time steps (hours)

    FUSION = inputs["FUSION"]
    dfFusion = inputs["dfFusion"]

    ##### Tritium Balance #####

    ## Tritium capacity
    @variable(EP, TritCap[y in FUSION] .>= 0.0)

    ## Tritium inventory
    @variable(EP, vtrit_inventory[y in FUSION, t=1:T], lower_bound = 0.)

    ## Tritium inventory <= Tritium capacity
    @constraint(EP, [y in FUSION,t=1:T], vtrit_inventory[y,t] <= TritCap[y])

    # Tritium exports
    @variable(EP, vtrit_exports[y in FUSION, t=1:T] >= 0.0)

    ## Hourly tritium decay rate [/ hr]
    TRIT_DECAY = 0.00000642259

    ## Fuel useage
    # Energy per fusion = 14e6 * 1.6e-19
    # Fusions per MWh = 3600 * 1e6 / (14e6 * 1.6e-19)
    # Avagadro's number = 6.022e23
    # Mol's per MWh = (3600 * 1e6 / (14e6 * 1.6e-19)) / 6.022e23
    # Tritium molar mass = 3g / mol
    # kg tritium per MWh = 3 / 1000 * (3600 * 1e6 / (14e6 * 1.6e-19)) / 6.022e23
    # Secondary boost = 15%
    # Adjusted kg tritium per MWh = 3 / 1000 * (3600 * 1e6 / (14e6 * 1.6e-19)) / 6.022e23 / 1.15 = 6.96295e-6

    TRIT_FUEL_RATIO = 3 / 1000 * (3600 * 1e6 / (14e6 * 1.6e-19)) / 6.022e23 / 1.15
    DEU_FUEL_RATIO = 2 / 1000 * (3600 * 1e6 / (14e6 * 1.6e-19)) / 6.022e23 / 1.15

    ## Hourly tritium balance
    @constraint(EP, [y in FUSION, t=1:T], vtrit_inventory[y,t] == vtrit_inventory[y, cyclicindex(t-1, T)] * (1 - TRIT_DECAY) + TRIT_FUEL_RATIO * (EP[:vThermOutput][y,t] * (dfFusion[y,:Trit_Breed] - dfFusion[y,:Trit_Fuel]) - EP[:num_units][y] * dfFusion[y,:Trit_Loss] - vtrit_exports[y,t]))

    ##### Deuterium Balance #####
    @variable(EP, DeuCap[y in FUSION] .>= 0.0)

    # Deuterium inventory
    @variable(EP, vdeu_inventory[y in FUSION, t=1:T], lower_bound = 0.)
    @constraint(EP, [y in FUSION,t=1:T], vdeu_inventory[y,t] <= DeuCap[y])

    # Deuterium exports
    @variable(EP, vdeu_exports[y in FUSION, t=1:T] >= 0.0)

    # Deuterium balance
    @constraint(EP, [y in FUSION, t=1:T], vdeu_inventory[y,t] == vdeu_inventory[y, cyclicindex(t-1, T)] - DEU_FUEL_RATIO * (EP[:vThermOutput][y,t] * dfFusion[y,:Deu_Fuel] + EP[:num_units][y] * dfFusion[y,:Deu_Loss] + vdeu_exports[y,t]))

    ### Add Tritium and Deuterium Storage Capacity to the Fixed Costs (Dummy Cost of $0.1/kg currently in place)

    # First, pull values of cost from dfFusion
    trit_stor_cost = dfFusion[!,:Trit_Stor_Cost]
    deu_stor_cost = dfFusion[!,:Deu_Stor_Cost]

    # Then, multiply by the storage capacity and add to the fixed costs
    EP[:eCFix][FUSION] .+= (TritCap[FUSION].*trit_stor_cost[FUSION].+DeuCap[FUSION].*deu_stor_cost[FUSION])

    @expression(EP, eSumFuelFix, (sum((TritCap[y].*trit_stor_cost[y].+DeuCap[y].*deu_stor_cost[y] for y in FUSION))))

    EP[:eTotalCFix] += eSumFuelFix

    EP[:eObj] += eSumFuelFix 

end

# This function incorporates the vessel costs
function fusionvessel(EP::Model, inputs::Dict, setup::Dict)
    dfGen = inputs["dfGen"]

    T = inputs["T"]     # Number of time steps (hours)

    FUSION = inputs["FUSION"]
    dfFusion = inputs["dfFusion"]

    # Grab variables from the fusion_data.csv
    plant_life = dfFusion[!,:Plant_Life]
    vessel_name = dfFusion[!,:Vessel_Life]
    replace_dur = dfFusion[!,:Rep_Dur]
    discount = dfFusion[!,:Dis_Rate]
    vessel_inv = dfFusion[!,:Inv_Vessel_per_MWe]
    plant_cost = dfGen[!,:Inv_Cost_per_MWyr]

    # Turbine efficiency
    MMBTU_PER_MWH = 3.412
    con_mWt = MMBTU_PER_MWH./dfGen[!,:Heat_Rate_MMBTU_per_MWh]

    ## Calculate annual plant utilization
    @expression(EP, eThermOutputTot[y in FUSION], sum(EP[:vThermOutput][y,t] for t=1:T))

    ## Calculates the effective thermal capacity of the plant for the year
    @expression(EP, eECap[y in FUSION], (8760 * dfGen[!,:Cap_Size][y]* (1/con_mWt[y]) * (20. / 21.)))

    ## Calculates Annual Utilization based on thermal output and capacity
    @expression(EP, eAnnualUtil[y in FUSION], eThermOutputTot[y]/eECap[y])

    ## Constrain annual utilization
    @constraint(EP, cAnnualUtil[y in FUSION], eAnnualUtil[y] <= ((-vessel_name[y] + sqrt((vessel_name[y])^2 + 4*vessel_name[y]*replace_dur[y])) / (2*replace_dur[y])))

    ## Annuitized Plant costs
    @expression(EP, ePlantAnnual[y in FUSION], (plant_cost[y] * discount[y])/(1 - (1 + discount[y])^(-plant_life[y])))

    ## Calculate Vessel Investment Costs
    @expression(EP, eC1[y in FUSION], vessel_inv[y] * ((0.5*discount[y] / ((1+discount[y])^(2*vessel_name[y]) - 1)) - ((2*vessel_name[y]*discount[y]*(1 + discount[y])^(2*vessel_name[y]) * log(1+discount[y])) / ((1+discount[y])^(2*vessel_name[y]) -1)^2)))

    @expression(EP, eC2[y in FUSION], vessel_inv[y] * ((4*vessel_name[y]*discount[y]*(1+discount[y])^(2*vessel_name[y]) * log(1+discount[y])) / ((1+discount[y])^(2*vessel_name[y]) - 1)^2))

    ## Fixed Vessel Investment Costs
    @expression(EP, eVesselFix[y in FUSION], (vessel_inv[y]*discount[y])/(1-(1+discount[y])^(-plant_life[y])))
    
    ## Add Vessel Investment Costs to Fixed/Var Costs
    EP[:eCFix][FUSION] .+= (eC1[FUSION].+eVesselFix[FUSION]).*EP[:eTotalCap][FUSION]
    EP[:eCVar_out][FUSION] .+= eC2[FUSION].*eThermOutputTot[FUSION].*(20. / 21.).*(1/8760.).*con_mWt[FUSION]

    ## Sum Vessel Investment and Variable Costs to add to TotalFix and TotalVar
    @expression(EP, eSumCFix, sum((eC1[y] + eVesselFix[y]) * EP[:eTotalCap][y] for y in FUSION))

    pulse_ratio = fusionpulseratio.(dfFusion[!,:Pulse_length_min], dfFusion[!,:Pulse_gap_min])
    @expression(EP, eSumCVar, sum(eC2[y] * eThermOutputTot[y] * con_mWt[y] * pulse_ratio[y] for y in FUSION) / T)

    EP[:eTotalCFix] += eSumCFix
    EP[:eTotalCVarOut] += eSumCVar

    ## Add the values to eObj
    EP[:eObj] += eSumCFix + eSumCVar
end

function fusionthermalstorage(EP::Model, inputs::Dict, setup::Dict)
    T = inputs["T"]     # Number of time steps (hours)

    FUSION = inputs["FUSION"]
    dfFusion = inputs["dfFusion"]

    ## Grab cost inputs from dfFusion
    stor_cost = dfFusion[!,:Stor_Cost_per_MWht]
    dis_cost = dfFusion[!,:Dis_Cost_per_MWht]

    ## Thermal storage energy capacity [MWht]
    @variable(EP, vThermStorCap[y in FUSION], lower_bound = 0.)

    ## Thermal storage hourly charge and discharge capacity [MWht]
    @variable(EP, vThermDisCap[y in FUSION], lower_bound = 0.)
    
    ## Discharge per hour [MWht]
    @variable(EP, vThermDis[y in FUSION, t=1:T], lower_bound = 0.)
    
    ## Charge per hour [MWht]
    @variable(EP, vThermChar[y in FUSION, t=1:T], lower_bound = 0.)

    ## Energy stored in thermal storage [MWht]
    @variable(EP, vThermStor[y in FUSION, t=1:T], lower_bound = 0.)

    ## Energy stored <= Storage Capacity
    @constraint(EP, [y in FUSION,t=1:T], vThermStor[y,t] <= vThermStorCap[y])

    ## Limit charge and discharge rates
    @constraint(EP, [y in FUSION,t=1:T], vThermDis[y,t]  <= vThermDisCap[y])
    @constraint(EP, [y in FUSION,t=1:T], vThermChar[y,t] <= vThermDisCap[y])

    ## Calculate the net discharge rate
    @expression(EP, eThermStorNetDischarge[y in FUSION, t=1:T], vThermDis[y,t] - vThermChar[y,t])

    ## Change in stored energy == net discharge rate
    ## Later, we can add losses or other factors which differentiate these two terms
    @constraint(EP, [y in FUSION,t=1:T], vThermStor[y,cyclicindex(t-1,T)] - vThermStor[y,t] == eThermStorNetDischarge[y,t])

    ## Fixed costs associated with the thermal storage
    thermstoragecosts = vThermStorCap[FUSION] .* stor_cost[FUSION] .+ vThermDisCap[FUSION] .* dis_cost[FUSION]

    ## Add thermal storage fixed costs to the fusion generator costs
    EP[:eCFix][FUSION] .+= thermstoragecosts

    @expression(EP, eSumThermFix, sum(thermstoragecosts[y] for y in FUSION))

    ## Add the fusion fleets thermal storage fixed costs to total grid costs
    EP[:eTotalCFix] += eSumThermFix

    ## Add the fusion fleets thermal storage fixed costs to the objective function
    EP[:eObj] += eSumThermFix
end

function cyclicindex(idx::Int, maxidx::Int)
    return mod(idx - 1, maxidx) + 1
end

function fusionpulseratio(pulse_length_min::Float64, pulse_gap_min::Float64)
    return pulse_length_min / (pulse_length_min + pulse_gap_min)
end

## This function is the overall function for fusion power 
function fusion!(EP::Model, inputs::Dict, setup::Dict)
    fusionthermalpower(EP,inputs,setup)
    fusionthermalstorage(EP,inputs,setup)
    fusiongridpower(EP,inputs,setup)
    fusionfuel(EP,inputs,setup)
    fusionvessel(EP,inputs,setup)
end     