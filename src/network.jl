struct Network
    nPopulations::Int
    k_bar::Int
end

function fillConnectionMatrix!(connections::Array{Float64, 2}, nPopulations::Int, k_bar)
    for i in 1:nPopulations
        columns = [mod(i+j-1, nPopulations)+1 for j in -k_bar÷2:k_bar÷2 if mod(i+j-1, nPopulations)+1 != i]
        connections[i, columns] .= 1
    end
end


function updateNetwork!(populations, connections,epi)
    populations_copy = deepcopy(populations)

    for (i, population_copy) in enumerate(populations_copy)
        updatePopulation!(populations[i], connections[i,:],populations_copy,epi)
        # populations[i] = population_copy
    end
end

function computePOConnectivityHistory(mobilityRatesHistory,infectedHistory,populations,connections,POrder)
    # POrder =0, 1, or 2
    nTimeSteps, nPopulations, nMobility = size(mobilityRatesHistory)
    POConnectivity = zeros(nTimeSteps, nPopulations)
    AverageP2Connectivity = zeros(nTimeSteps)
    for t in 1:nTimeSteps
        for i in 1:nPopulations
            for j in 1:nMobility
                POConnectivity[t, i] += mobilityRatesHistory[t, i, j]*mobilityRatesHistory[t, j, i]*connections[i,j]
                if POrder > 0 POConnectivity[t, i] *= (populations[j].size - infectedHistory[t,i]) end 
                if POrder > 1 POConnectivity[t, i] *= (populations[i].size-infectedHistory[t,j]) end
            end
            AverageP2Connectivity[t]+=POConnectivity[t, i]
        end
    end
    AverageP2Connectivity = AverageP2Connectivity./nPopulations
    return POConnectivity, AverageP2Connectivity
end
