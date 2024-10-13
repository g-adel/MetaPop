mutable struct Population
    S::Float64 # Susciptible population fraction
    I::Float64 # Infected population fraction
    R::Float64 # Recovered population fraction
    size::Float64 # Total population size
    position::Point # Position of the population
    index::Int # Index number of the population
    ρs::Array{Float64, 1} # Mobility rates to other populations
end

@kwdef mutable struct Metapopulation
    S
    populations::Array{Population, 1}
    mobilityRates::Array{Float64, 2} 
end

mutable struct PopulationRoC
    dS::Float64 # Susciptible population fraction RoC
    dI::Float64 # Infected population fraction RoC
    dR::Float64 # Recovered population fraction RoC
    ρsRoC::Array{Float64, 1}
end

function initializePopulations!(meta)
    r = W/2 - 30
    populations=meta.populations
    nPopulations = meta.S.net.nPopulations
    for i in 1:nPopulations
        theta = 2*π/nPopulations * (i-1)
        location = Point(r*cos(theta), r*sin(theta))
        populations[i] = Population(1., 0., 0., 1., location, i,zeros(nPopulations))
    end
    populations[1].I=meta.S.sim.I₀;
    populations[1].S= 1 - populations[1].I
    
end


function getPopulationRoC(pop::Population,meta::Metapopulation)
    net, populations, epi = meta.S.net, meta.populations, meta.S.epi;
    S = pop.S; I = pop.I; R = pop.R
    inConnections = net.connections[:,pop.index]
    outConnections = net.connections[pop.index,:]
    netFlowInfected = netFlowSusceptible = netFlowRecovered = 0
    for (connPopInd, connWeight) in enumerate(inConnections)
        connPop = populations[connPopInd]
        finalMobilityRate = meta.mobilityRates[pop.index,connPopInd]

        netFlowSusceptible += finalMobilityRate * (connPop.S)
        netFlowInfected += finalMobilityRate * (connPop.I)
        netFlowRecovered += finalMobilityRate * (connPop.R)
    end
    for (connPopInd, connWeight) in enumerate(outConnections)
        connPop = populations[connPopInd]
        finalMobilityRate = meta.mobilityRates[connPopInd, pop.index]

        netFlowSusceptible += finalMobilityRate * (- S)
        netFlowInfected += finalMobilityRate    * (- I)
        netFlowRecovered += finalMobilityRate   * (- R)
    end
    dS = -epi.β*I*S + epi.σ*R + netFlowSusceptible # is netFlow a rate of change?
    dI =  epi.β*I*S - epi.γ*I + netFlowInfected
    dR =  epi.γ*I   - epi.σ*R + netFlowRecovered
    strategy = meta.S.strat.strategy
    ρsRoC = empty(size(populations))
    if strategy == GlobalDiffRestriction
        ρsRoC == globalDiffRestriction(pop,inConnections,meta)
    elseif strategy == UniformDiffRestriction
        ρsRoC = uniformDiffRestriction(pop,inConnections,meta)
    elseif strategy == IndivDiffRestriction
        ρsRoC = indivDiffRestriction(pop,inConnections,meta)
    end
    populationRoC = PopulationRoC(dS,dI,dR,ρsRoC) #struct
    return populationRoC
end