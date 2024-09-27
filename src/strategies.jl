Base.@kwdef struct Strat
    λ::Float64
    mobBias::Float64
end

function uniformDiffRestriction(pop,localConnections::Array{Float64, 1},meta)
    λ = pop.strat.λ # Adaptive mobility tuning rate
    nbrs_indices = findall(x -> x > 0, localConnections) # TODO: store as sparse matrix
    
    netFlowInfected = 0;
    for connPopInd in nbrs_indices
        # could be turned into array operation
        finalMobilityRate = meta.mobilityRates[connPopInd,pop.index]
        netFlowInfected += finalMobilityRate * (meta.populations[connPopInd].I)
    end

    restrictionRoC = zeros(size(pop.restrictions))
    for i in nbrs_indices
        restrictionRoC[i] = λ*netFlowInfected - pop.strat.mobBias * pop.restrictions[i]
    end

    return restrictionRoC
end