struct Strat
    λ::Float64
    mobBias::Float64
end

# return array of floats of size nPopulations for mobility rate 

function uniformDiffRestriction(pop,populations,localConnections::Array{Float64, 1},epi)
    # Return new mobility rates
    λ = pop.strat.λ # Adaptive mobility tuning rate
    _, _, _, μ = structVals(epi)
    # compute average incoming and outgoing infection rate
    netFlowInfected = 0;
    nbrs_indices = findall(x -> x > 0, localConnections) # TODO: store as sparse matrix
    
    for (connPopInd, connWeight) in enumerate(localConnections)        
        connPop = populations[connPopInd]
        finalMobilityRate = connWeight * μ *  (1-connPop.mobilityRestrictions[pop.index]) * (1-pop.mobilityRestrictions[connPopInd])

        netFlowInfected += finalMobilityRate * (populations[connPopInd].I)# -pop.I)
        # println("connPopInd: ", connPopInd, ", connWeight: ", connWeight, ", finalMobilityRate: ", finalMobilityRate, ", netFlowInfected: ", netFlowInfected)
    end

    obilityRestrictionResidual = λ*netFlowInfected - pop.strat.mobBias ### STRATEGY
    
    newMobilityRestrictions = deepcopy(pop.mobilityRestrictions)    
    for i in nbrs_indices
        newMobilityRestrictions[i] = pop.mobilityRestrictions[i] + obilityRestrictionResidual
    end
    clamp!(newMobilityRestrictions,0.0,.99)
    return newMobilityRestrictions
end

# function homogenousLinearRestriction(pop,populations,localConnections::Array{Float64, 1},epi)

### NEEDS TO BE REWRITTEN
# function linearIndivDiffRestriction(popInd,populations,localConnections::Array{Float64, 1},epi; bias=0)
#     # Return new mobility rates
#     λ = 8000 # Adaptive mobility tuning rate
#     γ, β, μ = structVals(epi)
#     pop = populations[popInd]
#     # compute average incoming and outgoing infection rate
#     newMobilityRestrictions = deepcopy(pop.mobilityRestrictions)
#     mobilityBias = bias/λ*(pop.size-pop.I)
#     for (connPopInd, connWeight) in enumerate(localConnections)
#         connPop = populations[connPopInd]
#         if connWeight > 0
#             linkMobilityRate = μ * connWeight * connPop.mobilityRestrictions[popInd] * pop.mobilityRestrictions[connPopInd]
#             incomingInfected = linkMobilityRate * populations[connPopInd].I
#             outgoingInfected = linkMobilityRate * pop.I
#             # (connPopInd == 10 && popInd == 8) &&   print("incomingInfected: ",incomingInfected,", outgoingInfected: ",outgoingInfected,"\n")
#             newMobilityRestrictions[connPopInd] += λ*(outgoingInfected  - incomingInfected + mobilityBias)*pop.mobilityRestrictions[connPopInd]
#         end
#     end
#     clamp!(newMobilityRestrictions,0.001,1.)
#     return newMobilityRestrictions
# end