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

function write_start(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
	dfGen = inputs["dfGen"]
	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)
	COMMIT = inputs["COMMIT"]
	# Startup state for each resource in each time step
	dfStart = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!, :Zone])
	start = zeros(T,G)
	start[:, COMMIT] = value.(EP[:vSTART][:, COMMIT])
	dfStart.AnnualSum = start' * inputs["omega"]
	dfStart = hcat(dfStart, DataFrame(start', :auto))
	auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("AnnualSum");[Symbol("t$t") for t in 1:T]]
	rename!(dfStart,auxNew_Names)

	total = DataFrame(["Total" 0 sum(dfStart.AnnualSum) fill(0.0, (1,T))], :auto)
	total[:, 4:T+3] .= sum(start', dims = 1)
	rename!(total,auxNew_Names)
	dfStart = vcat(dfStart, total)
	CSV.write(joinpath(path, "start.csv"), dftranspose(dfStart, false), writeheader=false)
end

function write_start_3(path::AbstractString, inputs::Dict{String, Any}, setup::Dict{Any, Any}, EP::Model)
	dfGen = inputs["dfGen"]::DataFrame
	G = inputs["G"]::Int64   # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]::Int64   # Number of time steps (hours)
	COMMIT = inputs["COMMIT"]::Vector{Int64}
	# Startup state for each resource in each time step

	resources = inputs["RESOURCES"]::Vector{String31}
	dfStart = DataFrame([[] for _ = resources], resources)::DataFrame
	push!(dfStart, dfGen[!, :Zone])
	
	start = zeros(T,G)
	start[:, COMMIT] = value.(EP[:vSTART][:, COMMIT])
	push!(dfStart, inputs["omega"]' * start)

	dfStart = vcat(dfStart, DataFrame(start, resources))

	insertcols!(dfStart, "Total"  => sum.(eachrow(dfStart)))
	dfStart.Total[1] = 0

	auxNew_Names=[Symbol("Zone");Symbol("AnnualSum");[Symbol("t$t") for t in 1:T]]
	insertcols!(dfStart, 1, "Resource"  => auxNew_Names)
	
	CSV.write(joinpath(path, "start.csv"), dfStart, writeheader=true)
end