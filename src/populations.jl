mutable struct Population
    S::Float64
    I::Float64
    R::Float64
    size::Float64
    position::Point
    index::Int
    mobilityRates::Array{Float64, 1}
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
    S = pop.S; I = pop.I; R = pop.R
    α, β, μ = structVals(epi)
    totalIncomingInfected = totalIncomingSusceptible = totalIncomingRecovered = 0
    totalOutgoingInfected = totalOutgoingSusceptible = totalOutgoingRecovered = 0
    for (connPopInd, connectionWeight) in enumerate(localConnections)
        if connectionWeight > 0
            finalMobilityRate = μ * connectionWeight * populations[connPopInd].mobilityRates[pop.index] * populations[pop.index].mobilityRates[connPopInd]
            totalIncomingSusceptible += finalMobilityRate * populations[connPopInd].S
            totalIncomingInfected += finalMobilityRate * populations[connPopInd].I
            totalIncomingRecovered += finalMobilityRate * populations[connPopInd].R
            
            totalOutgoingSusceptible += finalMobilityRate * S
            totalOutgoingInfected += finalMobilityRate * I
            totalOutgoingRecovered += finalMobilityRate * R
        end
    end
    dS = -(β*I)*S + totalIncomingSusceptible - totalOutgoingSusceptible
    dI = (β*I)*S - α*I + totalIncomingInfected - totalOutgoingInfected
    dR = α*I + totalIncomingRecovered - totalOutgoingRecovered

    # integration
    pop.S += dS
    pop.I += dI
    pop.R += dR

    # Mobility adaptivity update
    # pop.mobilityRates = linearUniformDiffResponse(pop.index, populations,localConnections,epi,bias = 0)
    pop.mobilityRates = linearIndivDiffResponse(pop.index, populations,localConnections,epi,bias = 0)
end
