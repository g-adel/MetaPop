@enum StrategyType begin
    GlobalDiffRestriction
    UniformDiffRestriction
    IndivDiffRestriction
    IndivLogRestriction
end

Base.@kwdef struct Strat
    λ::Float64
    mobBias::Float64
    strategy::StrategyType
end

function globalDiffRestriction(pop,localConnections::Array{Float64,1},meta)
end

function uniformPropRestriction(pop,localConnections::Array{Float64, 1},meta)
    λ = meta.S.strat.λ # Adaptive mobility tuning rate
    nbrs_indices = findall(x -> x > 0, localConnections) # TODO: store as sparse matrix
    
    netFlowInfected = 0;
    for connPopInd in nbrs_indices
        # could be turned into array operation
        finalMobilityRate = meta.mobilityRates[connPopInd,pop.index]
        netFlowInfected += finalMobilityRate * (meta.populations[connPopInd].I)
    end

    ρRoC = zeros(size(pop.ρs))
    for i in nbrs_indices
        ρRoC[i] = λ*netFlowInfected - meta.S.strat.mobBias * pop.ρs[i]
    end

    return ρRoC
end

function indivPropRestriction(pop,localConnections::Array{Float64, 1},meta)
    λ = meta.S.strat.λ # Adaptive mobility tuning rate
    nbrs_indices = findall(x -> x > 0, localConnections) # TODO: store as sparse matrix
    ρRoC = zeros(size(pop.ρs))
    
    for connPopInd in nbrs_indices
        # could be turned into array operation
        inflowInfected = meta.mobilityRates[connPopInd,pop.index] * (meta.populations[connPopInd].I)
        ρRoC[connPopInd] = λ*inflowInfected - meta.S.strat.mobBias * pop.ρs[connPopInd]
    end

    return ρRoC
end

function indivLogRestriction(pop,localConnections::Array{Float64, 1},meta)
    λ = meta.S.strat.λ # Adaptive mobility tuning rate
    nbrs_indices = findall(x -> x > 0, localConnections) # TODO: store as sparse matrix
    ρRoC = zeros(size(pop.ρs))
    
    for connPopInd in nbrs_indices
        # could be turned into array operation
        inflowInfected = meta.mobilityRates[connPopInd,pop.index] * (meta.populations[connPopInd].I)
        localRate = (meta.S.epi.β - meta.S.epi.γ) * pop.I
        # ρRoC[connPopInd] = λ*log((inflowInfected + localRate)/localRate) - meta.S.strat.mobBias * pop.ρs[connPopInd]
        ρRoC[connPopInd] = λ*log(meta.populations[connPopInd].I/pop.I) - meta.S.strat.mobBias * pop.ρs[connPopInd]

        if isnan(ρRoC[connPopInd]) || isinf(ρRoC[connPopInd])
            ρRoC[connPopInd] = 0
        end
    end

    return ρRoC
end