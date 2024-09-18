mutable struct Population
    S::Float64 # Susciptible population fraction
    I::Float64 # Infected population fraction
    R::Float64 # Recovered population fraction
    size::Float64 # Total population size
    position::Point # Position of the population
    index::Int # Index number of the population
    strat#::Strat # strategy parameters
    restrictions::Array{Float64, 1} # Mobility rates to other populations
end

mutable struct PopulationRoC
    dS::Float64 # Susciptible population fraction
    dI::Float64 # Infected population fraction
    dR::Float64 # Recovered population fraction
    restrictionsRoC::Array{Float64, 1}
end

function initializePopulations!(populations,strat)
    r = W/2 - 30
    nPopulations = length(populations)
    for i in 1:nPopulations
        theta = 2*π/nPopulations * (i-1)
        location = Point(r*cos(theta), r*sin(theta))
        populations[i] = Population(1., 0., 0., 1., location, i,strat,zeros(nPopulations))
    end
    populations[1].I=.00003;
    populations[1].S= 1 - populations[1].I
    
end


function getPopulationRoC(pop::Population,localConnections::Array{Float64, 1},populations,epi)
    S = pop.S; I = pop.I; R = pop.R
    γ, β, σ, μ = structVals(epi)
    netFlowInfected = netFlowSusceptible = netFlowRecovered = 0
    for (connPopInd, connWeight) in enumerate(localConnections)
        connPop = populations[connPopInd]
        finalMobilityRate = connWeight * μ *  (1-connPop.restrictions[pop.index]) * (1-connPop.restrictions[connPopInd])

        netFlowSusceptible += finalMobilityRate * (connPop.S - S)
        netFlowInfected += finalMobilityRate * (connPop.I - I)
        netFlowRecovered += finalMobilityRate * (R - connPop.R)
    end
    dS = -β*I*S + σ*R + netFlowSusceptible # is netFlow a rate of change?
    dI = β*I*S  - γ*I + netFlowInfected
    dR = γ*I    - σ*R + netFlowRecovered

    # Restriction RoC
    restrictionsRoC = uniformDiffRestriction(pop, populations,localConnections,epi)
    # restrictionsRoC = zeros(size(populations))
    populationRoC = PopulationRoC(dS,dI,dR,restrictionsRoC)
    return populationRoC
end