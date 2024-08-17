

# return array of floats of size nPopulations for mobility rate 

function linearUniformDiffResponse(popInd,populations,localConnections::Array{Float64, 1},epi; bias=0)
    # Return new mobility rates
    λ = 8000 # Adaptive mobility tuning rate
    γ, β, μ = structVals(epi)
    pop = populations[popInd]
    # compute average incoming and outgoing infection rate
    totalIncomingInfected = 0;    totalOutgoingInfected = 0; 
    totConnWeight = 0
    
    for (connPopInd, connWeight) in enumerate(localConnections)
        finalMobilityRate = connWeight * μ  * populations[connPopInd].mobilityRates[popInd] * populations[popInd].mobilityRates[popInd]
        totalIncomingInfected += finalMobilityRate * populations[connPopInd].I
        totalOutgoingInfected += finalMobilityRate * pop.I
        totConnWeight+=connWeight
    end
    netRateAvg = (totalOutgoingInfected - totalIncomingInfected)/totConnWeight
    newMobilityRates = deepcopy(pop.mobilityRates)    
    mobilityBias = bias/λ *(pop.size-pop.I)

    for (i, mobilityRate) in enumerate(newMobilityRates)
        newMobilityRates[i] += (localConnections[i]>0)*(netRateAvg + mobilityBias)*pop.mobilityRates[i]*λ
    end
    clamp!(newMobilityRates,0.001,1.)
    # println("Local Connection I ",localConnections[popInd]," newMobilityRates I ",newMobilityRates[popInd])
    return newMobilityRates
end

function linearIndivDiffResponse(popInd,populations,localConnections::Array{Float64, 1},epi; bias=0)
    # Return new mobility rates
    λ = 8000 # Adaptive mobility tuning rate
    γ, β, μ = structVals(epi)
    pop = populations[popInd]
    # compute average incoming and outgoing infection rate
    totConnWeight = 0
    newMobilityRates = deepcopy(pop.mobilityRates)
    mobilityBias = bias/λ*(pop.size-pop.I)
    for (connPopInd, connWeight) in enumerate(localConnections)
        if connWeight > 0
            totConnWeight+=connWeight
            linkMobilityRate = μ * connWeight * populations[connPopInd].mobilityRates[popInd] * populations[popInd].mobilityRates[connPopInd]
            incomingInfected = linkMobilityRate * populations[connPopInd].I
            outgoingInfected = linkMobilityRate * pop.I
            # (connPopInd == 10 && popInd == 8) &&   print("incomingInfected: ",incomingInfected,", outgoingInfected: ",outgoingInfected,"\n")
            newMobilityRates[connPopInd] += λ*(outgoingInfected  - incomingInfected + mobilityBias)*pop.mobilityRates[connPopInd]
        end
    end
    clamp!(newMobilityRates,0.001,1.)
    return newMobilityRates
end