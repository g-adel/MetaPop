mutable struct Population
    S::Float64 # Susciptible population fraction
    I::Float64 # Infected population fraction
    R::Float64 # Recovered population fraction
    size::Float64 # Total population size
    position::Point # Position of the population
    index::Int # Index number of the population
    mobilityRates::Array{Float64, 1} # Mobility rates to other populations
end

function initializePopulations!(populations)
    r = W/2 - 30
    nPopulations = length(populations)
    uniformMobilityRates = ones(nPopulations)
    for i in 1:nPopulations
        theta = 2*π/nPopulations * (i-1)
        populations[i] = Population(1., 0., 0., 1., Point(r*cos(theta), r*sin(theta)), i,deepcopy(uniformMobilityRates))
    end
    populations[1].I=.0001;
    populations[1].S= 1 - populations[1].I
    
end


function updatePopulation!(pop::Population,localConnections::Array{Float64, 1},populations,epi)
    # DO NOT use populations[pop.index] equivalently with pop - not the same for order of updates
    S = pop.S; I = pop.I; R = pop.R
    γ, β, μ = structVals(epi)
    netFlowInfected = netFlowSusceptible = netFlowRecovered = 0
    for (connPopInd, connectionWeight) in enumerate(localConnections)
        finalMobilityRate = connectionWeight * μ *  populations[connPopInd].mobilityRates[pop.index] * populations[pop.index].mobilityRates[connPopInd]

        netFlowSusceptible += finalMobilityRate * (populations[connPopInd].S - S)
        netFlowInfected += finalMobilityRate * (populations[connPopInd].I - I)
        netFlowRecovered += finalMobilityRate * (R - populations[connPopInd].R)
    end
    dS = -(β*I)*S +  netFlowSusceptible
    dI = (β*I)*S - γ*I + netFlowInfected
    dR = γ*I + netFlowRecovered

    # integration
    pop.S += dS
    pop.I += dI
    pop.R += dR

    # Mobility adaptivity update
    pop.mobilityRates = linearUniformDiffResponse(pop.index, populations,localConnections,epi,bias = 0.1)
    # pop.mobilityRates = linearIndivDiffResponse(pop.index, populations,localConnections,epi,bias = 0.0)
end
