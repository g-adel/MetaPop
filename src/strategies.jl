Base.@kwdef struct Strat
    λ::Float64
    mobBias::Float64
end

# return array of floats of size nPopulations for mobility rate 

function globalDiffRestriction(pop,populations,globalInfectedFlow,epi)
    # Return restriction RoC
    λ = pop.strat.λ # Adaptive mobility tuning rate
    # compute average incoming and outgoing infection rate
    netFlowInfected = 0;
    nbrs_indices = findall(x -> x > 0, localConnections)
    
    for (connPopInd, connWeight) in enumerate(localConnections)
        connPop = populations[connPopInd]
        finalMobilityRate = connWeight * epi.μ *  max(0,(1-connPop.restrictions[pop.index])) * max(0,(1-pop.restrictions[connPopInd]))

        netFlowInfected += finalMobilityRate * (populations[connPopInd].I)
        # println("connPopInd: ", connPopInd, ", connWeight: ", connWeight, ", finalMobilityRate: ", finalMobilityRate, ", netFlowInfected: ", netFlowInfected)
    end

    # restrictionResidual = λ*netFlowInfected - pop.strat.mobBias *  ### STRATEGY
    restrictionRateOfChange = zeros(size(pop.restrictions))
    for i in nbrs_indices
        restrictionRateOfChange[i] = λ*netFlowInfected - pop.strat.mobBias * pop.restrictions[i]
    end
    return restrictionRateOfChange
end

function uniformDiffRestriction(pop,populations,localConnections::Array{Float64, 1},epi)
    # Return new mobility rates
    λ = pop.strat.λ # Adaptive mobility tuning rate
    _, _, _, μ = structVals(epi)
    # compute average incoming and outgoing infection rate
    netFlowInfected = 0;
    nbrs_indices = findall(x -> x > 0, localConnections) # TODO: store as sparse matrix
    
    for (connPopInd, connWeight) in enumerate(localConnections)        
        connPop = populations[connPopInd]
        finalMobilityRate = connWeight * μ *  max(0,(1-connPop.restrictions[pop.index])) * max(0,(1-pop.restrictions[connPopInd]))

        netFlowInfected += finalMobilityRate * (populations[connPopInd].I)
        # println("connPopInd: ", connPopInd, ", connWeight: ", connWeight, ", finalMobilityRate: ", finalMobilityRate, ", netFlowInfected: ", netFlowInfected)
    end

    # restrictionResidual = λ*netFlowInfected - pop.strat.mobBias *  ### STRATEGY
    restrictionRateOfChange = zeros(size(pop.restrictions))
    for i in nbrs_indices
        restrictionRateOfChange[i] = λ*netFlowInfected - pop.strat.mobBias * pop.restrictions[i]
    end
    return restrictionRateOfChange
end

# function homogenousLinearRestriction(pop,populations,localConnections::Array{Float64, 1},epi)

### NEEDS TO BE REWRITTEN
# function linearIndivDiffRestriction(popInd,populations,localConnections::Array{Float64, 1},epi; bias=0)
#     # Return new mobility rates
#     λ = 8000 # Adaptive mobility tuning rate
#     γ, β, μ = structVals(epi)
#     pop = populations[popInd]
#     # compute average incoming and outgoing infection rate
#     newRestrictions = deepcopy(pop.mobilityRestrictions)
#     mobilityBias = bias/λ*(pop.size-pop.I)
#     for (connPopInd, connWeight) in enumerate(localConnections)
#         connPop = populations[connPopInd]
#         if connWeight > 0
#             linkMobilityRate = μ * connWeight * connPop.mobilityRestrictions[popInd] * pop.mobilityRestrictions[connPopInd]
#             incomingInfected = linkMobilityRate * populations[connPopInd].I
#             outgoingInfected = linkMobilityRate * pop.I
#             # (connPopInd == 10 && popInd == 8) &&   print("incomingInfected: ",incomingInfected,", outgoingInfected: ",outgoingInfected,"\n")
#             newRestrictions[connPopInd] += λ*(outgoingInfected  - incomingInfected + mobilityBias)*pop.mobilityRestrictions[connPopInd]
#         end
#     end
#     clamp!(newRestrictions,0.001,1.)
#     return newRestrictions
# end