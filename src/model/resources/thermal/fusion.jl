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

    con_mWh_heatrt = dfGen[!,:Heat_Rate_MMBTU_per_MWh]./hr_unit   # Heat Rate of the technology is always greater than the mmbtu to mwh conversion
	## Decision variables for unit commitment
	# commitment state variable

    # Expression to represent the commitment state of the fusion power generators
    # Note, this is the number of commited units in the cluster if using linearized unit commitment
    @expression(EP, eFusionCommit[y in FUSION, t=1:T], EP[:vCOMMIT][y,t])
    # @constraint(EP, [y in FUSION, t=20:30], EP[:vCOMMIT][y,t] == 0)
    # @constraint(EP, [y in FUSION, t=1:T], eFusionCommit[y,t] >= 0.1)

    # ## Initialize the variable for the reactor being on/standby (CAN probably find a way to call the inputs to get this value - most likely the COMMIT variable)
    # @variable(EP, fusioncommit[t=1:T], Bin)

    ## Fusion Power Nameplate Capacity (this is gross cap size)
    # NameplateCap = (dfGen[!,:Cap_Size]) .* con_mWh_heatrt   ## Convert gross capacity to thermal capacity (still keeping it in MW)

    ## Convert nameplate capacity to effective thermal power capacity
    # Effective capacity is a constant fraction of the nameplate capacity
    # Formula is as follows: Effective Cap = NameplateCap*(Reactor Pulse Time/(Reactor Pulse Time + Dwell Time))
    # Dwell Time is 1 minute and reactor pulse time is 20 minutes
    # @expression(EP, eEffectiveCap[y in FUSION], NameplateCap[y]*(20/21))

    ## Reactor Thermal Output
    @variable(EP, vThermOutput[y in FUSION,t=1:T] >= 0)

    ## Constrain the Thermal Output with the Effective Capacity
    # @constraint(EP, [y in FUSION,t=1:T], vThermOutput[y,t] <= eEffectiveCap[y] * eFusionCommit[y,t])
    @constraint(EP, [y in FUSION,t=1:T], vThermOutput[y,t] <= eFusionCommit[y,t] * dfGen[!,:Cap_Size][y] * con_mWh_heatrt[y] * (20. / 21.))

    ## Force fusion power to be a part of the grid
    # @constraint(EP, [y in FUSION], EP[:eTotalCap][y] >= dfGen[!,:Cap_Size][y])
end

### This function will calculate the power provided by the reactor to the grid
function fusiongridpower(EP::Model, inputs::Dict, setup::Dict)
    T = inputs["T"]     # Number of time steps (hours)
    Z = inputs["Z"]     # Number of zones

    dfGen = inputs["dfGen"]
    dfFusion = inputs["dfFusion"]

    FUSION = inputs["FUSION"]
    
    @expression(EP, num_units[y in FUSION], EP[:eTotalCap][y] / dfGen[!,:Cap_Size][y])

    hr_unit = 3.412 ### hard coded mmbtu to mwh conversion
    # con_mWh_heatrt = (dfGen[!,:Heat_Rate_MMBTU_per_MWh])./hr_unit   # Heat Rate of the technology is always greater than the mmbtu to mwh conversion

    ## Define key variables
    magcool = dfFusion[!,:Mag_Cool]  #Amount of power (MW) being delivered to cool the magnets

    ### Calculation for the thermal balance of the salt loop
    salteff = dfFusion[!,:Salt_Eff]      ## Salt electric heating efficiency (already converts MW to thermal)
    saltLosses = dfFusion[!,:Salt_Loss]  ## Hourly Thermal losses from salt loop (2 MW/hr) for each reactor

    @expression(EP, eplantfix[y in FUSION,t=1:T], 10 * EP[:eFusionCommit][y,t])  # Fixed power being delivered to the power plant (Number multiplied with the binary variable for start)
    @expression(EP, eplantvar[y in FUSION,t=1:T], (1/12) * EP[:vThermOutput][y,t])  # Variable power being delivered to the power plant (Function of the ThermalOutput of the power plant)
    
    @variable(EP, 0 <= vsaltpwr[y in FUSION,t=1:T])  # Variable of the electric power being delivered to heat the salt
    @constraint(EP, [y in FUSION,t=1:T], vsaltpwr[y,t] <= num_units[y] * saltLosses[y] / salteff[y])

    ## Combine to make recirculating power expression
    @expression(EP, eRecircpwr[y in FUSION,t=1:T], num_units[y] * magcool[y] + eplantfix[y,t] + vsaltpwr[y,t] + eplantvar[y,t])
    @constraint(EP, eRecircpwr .>= 0)

    ## Expression for the thermal energy entering the turbine 
    # @expression(EP, eTurbThermal[y in FUSION,t=1:T], vsaltpwr[y,t] * salteff[y] + EP[:vThermOutput][y,t] - num_units[y] * saltLosses[y])
    # @constraint(EP, eTurbThermal .>= 0)
    @variable(EP, vTurbThermal[y in FUSION,t=1:T] >= 0.0)
    @constraint(EP, [y in FUSION,t=1], vTurbThermal[y,t] == EP[:vThermOutput][y,t] - EP[:vThermStor][y,t] - num_units[y] * saltLosses[y] + EP[:vsaltpwr][y,t] * salteff[y])
    @constraint(EP, [y in FUSION,t=2:T], vTurbThermal[y,t] == EP[:vThermStor][y,t-1] - EP[:vThermStor][y,t] + EP[:vThermOutput][y,t] - num_units[y] * saltLosses[y] + EP[:vsaltpwr][y,t] * salteff[y])

    ## Define the turbine efficiency
    # turbeff = 0.40 #(Change out for heat rate in gen_data)

    ## Expression for the amount of electric power coming out from the turbine

    # First, define the variable turbine gross electric power capacity
    @variable(EP, vTurbElecCap[y in FUSION], lower_bound = 0.)

    # Then, define the expression for the amount of electric power exiting turbine
    @expression(EP, eTurbElec[y in FUSION,t=1:T], eTurbThermal[y,t] ./ (dfGen[y,:Heat_Rate_MMBTU_per_MWh] ./ hr_unit))


    # Place constraints on the amount of energy exiting the turbine
    @constraint(EP, [y in FUSION,t=1:T], eTurbElec[y,t] <= vTurbElecCap[y])
    @constraint(EP, eTurbElec[y,t] .>= 0)

    # Place constraint on turbine minimum power
    @constraint(EP, [y in FUSION,t=1:T], eTurbElec[y,t] <= vTurbElecCap[y,t] * EP[:eFusionCommit][y,t])
    @constraint(EP, [y in FUSION,t=1:T], eTurbElec[y,t] >= 0.4 * vTurbElecCap[y,t] * EP[:eFusionCommit][y,t])

    # Then, add the costs of the turbine to the fixed costs
    turb_cost = dfFusion[!,:Turb_Cap_Cost]

    # Adding costs of variable turbine capacity to fixed costs
    EP[:eCFix][FUSION] .+= (vTurbElecCap[FUSION].*turb_cost[FUSION])
    @expression(EP, eSumTurbFix, (sum((vTurbElecCap[y].*turb_cost[y] for y in FUSION))))
    EP[:eTotalCFix] += eSumTurbFix
    EP[:eObj] += eSumTurbFix 

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

    ## Constraints for the electric power node
    # @constraint(EP, [y in FUSION, t=1:T], eFusionNetElec[y,t] == EP[:vP][y,t])
    @constraint(EP, PwrNode[y in FUSION, t=1:T], eTurbElec[y,t] + vfusionimports[y,t] == EP[:vP][y,t] + EP[:eRecircpwr][y,t])
end


function fusionfuel(EP::Model, inputs::Dict, setup::Dict)
    dfGen = inputs["dfGen"]

    T = inputs["T"]     # Number of time steps (hours)

    FUSION = inputs["FUSION"]
    dfFusion = inputs["dfFusion"]

    ##### Tritium Balance #####

    trit_stor_cap = dfFusion[!,:Tritium_Cap]

    # Tritium inventory
    # @variable(EP, vtrit_inventory[y in FUSION, t=1:T], lower_bound = 0., upper_bound = dfFusion[!,:Tritium_Cap][y])
    @variable(EP, vtrit_inventory[y in FUSION, t=1:T], lower_bound = 0.)


    # Tritium exports
    @variable(EP, vtrit_exports[y in FUSION, t=1:T] >= 0.0)

    # Fuel Coefficient, Breeding Coefficient, Loss Rate
    # @expression(EP, etrit_fuel[y in FUSION,t=1:T], dfFusion[y,:Trit_Fuel])
    # @expression(EP, etrit_breed[y in FUSION,t=1:T], dfFusion[y,:Trit_Breed])
    # @expression(EP, etrit_loss[y in FUSION,t=1:T], dfFusion[y,:Trit_Loss]

    # Tritium Balance
    # @constraint(EP, [y in FUSION, t=1], vtrit_inventory[y,t] == dfFusion[!,:Tritium_Cap][y] / 2.0)
    @constraint(EP, [y in FUSION, t=1],   vtrit_inventory[y,t] == vtrit_inventory[y,T]   + EP[:vThermOutput][y,T] * (dfFusion[y,:Trit_Breed] - dfFusion[y,:Trit_Fuel]) - EP[:num_units][y] * dfFusion[y,:Trit_Loss] - vtrit_exports[y,T])
    @constraint(EP, [y in FUSION, t=2:T], vtrit_inventory[y,t] == vtrit_inventory[y,t-1] + EP[:vThermOutput][y,t] * (dfFusion[y,:Trit_Breed] - dfFusion[y,:Trit_Fuel]) - EP[:num_units][y] * dfFusion[y,:Trit_Loss] - vtrit_exports[y,t])

    ##### Deuterium Balance #####
    deu_stor_cap = dfFusion[!,:Deu_Cap]

    # Deuterium inventory
    # @variable(EP, vdeu_inventory[y in FUSION, t=1:T], lower_bound = 0., upper_bound = dfFusion[!,:Deu_Cap][y])
    @variable(EP, vdeu_inventory[y in FUSION, t=1:T], lower_bound = 0.)

    # Deuterium exports
    @variable(EP, vdeu_exports[y in FUSION, t=1:T] >= 0.0)

    # Deuterium balance
    # @constraint(EP, [y in FUSION, t=1], vdeu_inventory[y,t] == dfFusion[!,:Deu_Cap][y] / 2.0)
    @constraint(EP, [y in FUSION, t=1],   vdeu_inventory[y,t] == vdeu_inventory[y,T]   - EP[:vThermOutput][y,T] * dfFusion[y,:Deu_Fuel] - EP[:num_units][y] * dfFusion[y,:Deu_Loss] - vdeu_exports[y,T])
    @constraint(EP, [y in FUSION, t=2:T], vdeu_inventory[y,t] == vdeu_inventory[y,t-1] - EP[:vThermOutput][y,t] * dfFusion[y,:Deu_Fuel] - EP[:num_units][y] * dfFusion[y,:Deu_Loss] - vdeu_exports[y,t])


    ### Add Tritium and Deuterium Storage Capacity to the Fixed Costs (Dummy Cost of $0.1/kg currently in place)

    # First, pull values of cost from dfFusion
    trit_stor_cost = dfFusion[!,:Trit_Stor_Cost]
    deu_stor_cost = dfFusion[!,:Deu_Stor_Cost]

    # Then, multiply by the storage capacity and add to the fixed costs
    EP[:eCFix][FUSION] .+= (trit_stor_cap[FUSION].*trit_stor_cost[FUSION].+deu_stor_cap[FUSION].*deu_stor_cost[FUSION])

    @expression(EP, eSumFuelFix, (sum((trit_stor_cap[y].*trit_stor_cost[y].+deu_stor_cap[y].*deu_stor_cost[y] for y in FUSION))))

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
    hr_unit = 3.412
    con_mWt = hr_unit./dfGen[!,:Heat_Rate_MMBTU_per_MWh]

    ## Calculate annual plant utilization
    @expression(EP, eThermOutputTot[y in FUSION], sum(EP[:vThermOutput][y,t] for t=1:T))

    @expression(EP, eECap[y in FUSION], (8760 * dfGen[!,:Cap_Size][y]* (1/con_mWt[y]) * (20. / 21.)))

    @expression(EP, eAnnualUtil[y in FUSION], eThermOutputTot[y]/eECap[y])

    # annual_util = 0

    # for t in 1:T
    #     hrly_util = (EP[:vThermOutput][8,t]/(8760 * dfGen[!,:Cap_Size][8] * (1/con_mWt[8]) * (20. / 21.)))
    #     annual_util += hrly_util
    # end

    # for y in inputs["FUSION"]
    #     sum()

    # @expression(EP, eAnnualUtil, annual_util)

    ## Constrain annual utilization
    @constraint(EP, cAnnualUtil[y in FUSION], eAnnualUtil[y] <= ((-vessel_name[y] + sqrt((vessel_name[y])^2 + 4*vessel_name[y]*replace_dur[y])) / (2*replace_dur[y])))

    ## Annuitized Plant costs
    @expression(EP, ePlantAnnual[y in FUSION], (plant_cost[y] * discount[y])/(1 - (1 + discount[y])^(-plant_life[y])))
    
    ## Add Plant Investment Costs to Fixed Cost of Fusion
    # EP[:eCFix][FUSION] .+= ePlantAnnual[FUSION].*dfGen[!,:Cap_Size][FUSION]

    ## Actual Vessel Life
    # @expression(EP, eVesselLife[y in FUSION], vessel_name[y]/eAnnualUtil[y] + replace_dur[y])


    ## Calculate Vessel Investment Costs
    @expression(EP, eC1[y in FUSION], vessel_inv[y] * ((0.5*discount[y] / ((1+discount[y])^(2*vessel_name[y]) - 1)) - ((2*vessel_name[y]*discount[y]*(1 + discount[y])^(2*vessel_name[y]) * log(1+discount[y])) / ((1+discount[y])^(2*vessel_name[y]) -1)^2)))

    @expression(EP, eC2[y in FUSION], vessel_inv[y] * ((4*vessel_name[y]*discount[y]*(1+discount[y])^(2*vessel_name[y]) * log(1+discount[y])) / ((1+discount[y])^(2*vessel_name[y]) - 1)^2))

    ## Fixed Vessel Investment Costs
    @expression(EP, eVesselFix[y in FUSION], (vessel_inv[y]*discount[y])/(1-(1+discount[y])^(-plant_life[y])))
    
    ## Variable Costs 
    # for y in inputs["FUSION"]
    #     sum(EP[:vThermOutput][y,t=1:T])
    # end

        

    
    # @expression(EP, eThermOutputTot[y in FUSION], 
    # for y in inputs["FUSION"]
    #     sum(EP[:vThermOutput][y,t=1:T])
    # end
    # )
    
    ## Add Vessel Investment Costs to Fixed/Var Costs
    EP[:eCFix][FUSION] .+= (eC1[FUSION].+eVesselFix[FUSION]).*EP[:eTotalCap][FUSION]
    EP[:eCVar_out][FUSION] .+= eC2[FUSION].*eThermOutputTot[FUSION].*(20. / 21.).*(1/8760.).*con_mWt[FUSION]

    ## Sum Vessel Investment and Variable Costs to add to TotalFix and TotalVar
    @expression(EP, eSumCFix, (sum((eC1[y]+eVesselFix[y])*EP[:eTotalCap][y] for y in FUSION)))
    @expression(EP, eSumCVar, (sum(eC2[y]*eThermOutputTot[y]*con_mWt[y] for y in FUSION).*(20. / 21.).*(1/8760.)))

    EP[:eTotalCFix] += eSumCFix
    EP[:eTotalCVarOut] += eSumCVar

    ## Add the values to eObj
    EP[:eObj] += eSumCFix + eSumCVar

    # @expression(EP, eCVesselInv[y in FUSION], ((vessel_inv[y]*discount[y]) / (1 - (1 + discount[y])^(-plant_life[y]))) - eC1[y] - eC2[y]*vAnnualUtil[y])

end

function fusionthermalstorage(EP::Model, inputs::Dict, setup::Dict)
    dfGen = inputs["dfGen"]

    T = inputs["T"]     # Number of time steps (hours)

    FUSION = inputs["FUSION"]
    dfFusion = inputs["dfFusion"]

    ## Grab cost inputs from dfFusion
    stor_cost = dfFusion[!,:Stor_Cost_per_MWht]
    dis_cost = dfFusion[!,:Dis_Cost_per_MWht]

    ## Create variables for storage cap and discharge cap
    @variable(EP, vThermStorCap[y in FUSION], lower_bound = 0.)
    @variable(EP, vThermDisCap[y in FUSION], lower_bound = 0.)

    ## Create a variable for the thermal storage over time
    @variable(EP, vThermStor[y in FUSION, t=1:T], lower_bound = 0.)
    @constraint(EP, [y in FUSION,t=1:T], vThermStor[y,t] <= vThermStorCap[y])

    ## Constrain the change in storage over time by discharge cap
    @constraint(EP, [y in FUSION, t=2:T], vThermStor[y,t] - vThermStor[y,t-1] <= vThermDisCap[y])
    @constraint(EP, [y in FUSION, t=2:T], -vThermDisCap[y] <= vThermStor[y,t] - vThermStor[y,t-1])

    ## Add the costs of the thermal storage to the fixed costs
    EP[:eCFix][FUSION] .+= (vThermStorCap[FUSION].*stor_cost[FUSION].+vThermDisCap[FUSION].*dis_cost[FUSION])

    @expression(EP, eSumThermFix, (sum((vThermStorCap[y].*stor_cost[y].+vThermDisCap[y].*dis_cost[y] for y in FUSION))))

    EP[:eTotalCFix] += eSumThermFix

    EP[:eObj] += eSumThermFix 
end


## This function is the overall function for fusion power 
function fusion!(EP::Model, inputs::Dict, setup::Dict)
    fusionthermalpower(EP,inputs,setup)
    fusionthermalstorage(EP,inputs,setup)
    fusiongridpower(EP,inputs,setup)
    fusionfuel(EP,inputs,setup)
    fusionvessel(EP,inputs,setup)
end     