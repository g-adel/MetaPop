struct Network
    nPopulations::Int
    k_bar::Int
end

function KRegRingMatrix(net::Network)
    connections::Array{Float64, 2} = zeros(Float64, net.nPopulations, net.nPopulations)
    for i in 1:net.nPopulations
        columns = [mod(i+j-1, net.nPopulations)+1 for j in -net.k_bar÷2:net.k_bar÷2 if mod(i+j-1, net.nPopulations)+1 != i]
        connections[i, columns] .= 1
    end
    return connections
end

function KRegChainMatrix(net::Network)
    connections::Array{Float64, 2} = zeros(Float64, net.nPopulations, net.nPopulations)
    for i in 1:net.nPopulations
        for j in max(1, i - net.k_bar ÷ 2):min(net.nPopulations, i + net.k_bar ÷ 2)
            if i != j
                connections[i, j] = 1
            end
        end
    end
    return connections
end

function KRegRingMatrix(net::Network)
    connections::Array{Float64, 2} = zeros(Float64, net.nPopulations, net.nPopulations)
    for i in 1:net.nPopulations
        columns = [mod(i+j-1, net.nPopulations)+1 for j in -net.k_bar÷2:net.k_bar÷2 if mod(i+j-1, net.nPopulations)+1 != i]
        connections[i, columns] .= 1
    end
    return connections
end

function updateNetwork!(populations, infectedFlows, connections,epi)
    populations_copy = deepcopy(populations)

    for (i, population_copy) in enumerate(populations_copy)
        updatePopulation!(populations[i], connections[i,:], populations_copy, epi)
        # populations[i] = population_copy
    end
    for (popInd, pop) in enumerate(populations)
        localConnections = connections[popInd,:]
        for (connPopInd, connWeight) in enumerate(localConnections)
            connPop = populations[connPopInd]
            infectedFlows[popInd,connPopInd] =  connWeight * epi.μ *  (1-connPop.mobilityRestrictions[pop.index]) * (1-connPop.mobilityRestrictions[connPopInd])
        end
    end
end

function computePOConnectivityHistory(mobilityRestrictionsHistory,infectedHistory,populations,connections,POrder)
    # POrder = 0, 1, or 2
    nTimeSteps, nPopulations, nMobility = size(mobilityRestrictionsHistory)
    POConnectivity = zeros(nTimeSteps, nPopulations)
    AveragePOConnectivity = zeros(nTimeSteps)
    for t in 1:nTimeSteps
        for i in 1:nPopulations
            for j in 1:nMobility
                POConnectivity[t, i] += mobilityRestrictionsHistory[t, i, j]*mobilityRestrictionsHistory[t, j, i]*connections[i,j]
                if POrder > 0 POConnectivity[t, i] *= (populations[j].size - infectedHistory[t,i]) end 
                if POrder > 1 POConnectivity[t, i] *= (populations[i].size-infectedHistory[t,j]) end
            end
            AveragePOConnectivity[t]+=POConnectivity[t, i]
        end
    end
    AveragePOConnectivity = AveragePOConnectivity./nPopulations
    return POConnectivity, AveragePOConnectivity
end
