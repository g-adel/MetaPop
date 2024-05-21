mutable struct Population
    S::Float64
    I::Float64
    R::Float64
    size::Float64
    position::Point
    index::Int
end

function initializePopulations!(populations)
    r = W/2 - 30
    for i in 1:nPopulations
        theta = 2*π/nPopulations * i
        populations[i] = Population(1, 0, 0, 1, Point(r*cos(theta), r*sin(theta)), i)
    end
    populations[1].I=.01;
    populations[1].S= 1 - populations[1].I
    
end


function updatePopulation!(pop::Population,localConnections::Array{Float64, 1},populations)
    S = pop.S; I = pop.I; R = pop.R
    totalIncomingInfected = 0; totalIncomingSusceptible = 0; totalIncomingRecovered = 0
    totalOutgoingInfected = 0; totalOutgoingSusceptible = 0; totalOutgoingRecovered = 0
    for (index, connectionStrength) in enumerate(localConnections)
        if connectionStrength > 0
            totalIncomingSusceptible += connectionStrength * μ * populations[index].S
            totalIncomingInfected += connectionStrength * μ * populations[index].I
            totalIncomingRecovered += connectionStrength * μ * populations[index].R
            
            totalOutgoingSusceptible += connectionStrength * μ * S
            totalOutgoingInfected += connectionStrength * μ * I
            totalOutgoingRecovered += connectionStrength * μ * R
        
        end
    end
    dS = -(β*I)*S + totalIncomingSusceptible - totalOutgoingSusceptible
    dI = (β*I)*S - α*I + totalIncomingInfected - totalOutgoingInfected
    dR = α*I + totalIncomingRecovered - totalOutgoingRecovered

    # integration
    pop.S += dS
    pop.I += dI
    pop.R += dR
end
