using CSV, DataFrames, Statistics

### Get Load Data

# load_df = CSV.read("Load_data.csv", DataFrame)
# hourly_load = load_df[!,:Load_MW_z1]
# total_load_kwhe = sum(hourly_load)*1000
# print(total_load_kwhe)

# ## Read Capacity and Power results and calculate utilization of fusion

# cap7_df = CSV.read("Results_7/capacity.csv", DataFrame)
# power7_df = CSV.read("Results_7/power.csv", DataFrame)

# capV_df = CSV.read("Fusion_Vessel_Fin/capacity.csv", DataFrame)
# powerV_df = CSV.read("Fusion_Vessel_Fin/power.csv", DataFrame)

# old_fusion_cap = cap7_df[end-1,:EndCap]
# old_fusion_power =  power7_df[3:end,"fusion"]
# old_fusion_util = (old_fusion_power ./ onet_cor)*100

# new_fusion_cap = capV_df[end-1,:EndCap]
# new_fusion_power = powerV_df[3:end,"fusion"]
# new_fusion_util = (new_fusion_power ./ net_cor)*100

# util_comp = DataFrame(Time = 1:8760)
# util_comp[!,"Old_Util"] = old_fusion_util
# util_comp[!,"New_Util"] = new_fusion_util

# print(maximum(new_fusion_util))
# print(maximum(old_fusion_util))

# CSV.write("Utilization.csv", util_comp)
# get_pwr = findall(x -> x > 10, fusion_power)
# fusion_actual = fusion_power[get_pwr]

# cap_fac_avg = mean(fusion_actual./fusion_cap)*100

# ## Compare electricity prices and utilization

# nprice_df = CSV.read("Fusion_Vessel_Fin/prices.csv", DataFrame)
# price_act = nprice_df[!,"1"]

# price_util_df = DataFrame(Time = 1:8760)
# price_util_df[!,"Util"] = new_fusion_util
# price_util_df[!,"Prices"] = price_act

# CSV.write("PriceUtil.csv", price_util_df)


#### Calculate Capacity Factor 

## First, find the maximum net electric power that corresponds to the maximum gross electric power
# nfusion_df = CSV.read("Fusion_Vessel_Fin/fusion.csv", DataFrame)

# gross_elec = nfusion_df[!,"Gross Electric"]
# gmax_indx = findfirst(==(maximum(gross_elec)), gross_elec)
# net_cor = nfusion_df[gmax_indx,"Net Electric"]

# println(gross_max)
# println(gmax_indx)
# println(net_cor)

# ## Find the maximum net electric power for the old case
# ofusion_df = CSV.read("Results_7/fusion.csv", DataFrame)
# ogross_elec = ofusion_df[!,"Gross Electric"]
# ogmax_indx = findfirst(==(maximum(ogross_elec)), ogross_elec)
# maximum(ogross_elec)
# onet_cor = ofusion_df[ogmax_indx,"Net Electric"]

# omax_elec = onet_cor*8760
# oact_elec = power7_df[2,"fusion"]
# ocap_fac = (oact_elec/omax_elec)*100

# ## Now that you have the maximum net_electric power, multiply it by 8760 to get the total power
# max_elec = net_cor*8760
# act_elec = powerV_df[2,"fusion"]

# cap_fac_f = (act_elec/max_elec)*100


### Utilization for Vessel and No Vessel case

vscap_df = CSV.read("Validation Runs/TE_Validation_Runs/Results_20/capacity.csv", DataFrame)
vspower_df = CSV.read("Validation Runs/TE_Validation_Runs/Results_20/power.csv", DataFrame)

vsnocap_df = CSV.read("Validation Runs/VS__Validation_Runs/Results_20/capacity.csv", DataFrame)
vsnpower_df = CSV.read("Validation Runs/VS__Validation_Runs/Results_20/power.csv", DataFrame)


vs_cap = vscap_df[end-1,:EndCap]
vs_power =  vspower_df[3:end,"fusion"]
vsutil = (vs_power ./ (vs_cap*(200.8/276)))*100

novs_cap = vsnocap_df[end-1,:EndCap]
novs_power = vsnpower_df[3:end,"fusion"]
novsutil = (novs_power ./ (novs_cap*(200.8/276)))*100

util_c = DataFrame(Time = 1:8760)
util_c[!,"With Vessel"] = vsutil
util_c[!,"Without Vessel"] = novsutil

println(mean(vsutil))
println(mean(novsutil))

CSV.write("VesselUtilizationCompare.csv", util_c)
