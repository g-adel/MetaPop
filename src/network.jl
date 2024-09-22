Base.@kwdef mutable struct Network
    nPopulations::Int
    k_bar::Int
    connections::Array{Float64,2}
    graph
end

function KRegRingMatrix(net::Network)
    connections::Array{Float64, 2} = zeros(Float64, net.nPopulations, net.nPopulations)
    for i in 1:net.nPopulations
        columns = [mod(i+j-1, net.nPopulations)+1 for j in -net.k_bar÷2:net.k_bar÷2 if mod(i+j-1, net.nPopulations)+1 != i]
        connections[i, columns] .= 1
    end
    return connections
end

function smallWorldMatrix(net::Network; β=.5)
    is_graph_connected = false
    g = SimpleGraph()
    while !is_graph_connected
        g = watts_strogatz(net.nPopulations, net.k_bar,β)
        is_graph_connected = Graphs.is_connected(g)
    end
    adj = convert(Matrix{Float64}, adjacency_matrix(g))
    return adj, g
end

function baraAlbert(net::Network)
    g = barabasi_albert(net.nPopulations, net.k_bar)
    adj = convert(Matrix{Float64}, adjacency_matrix(g))
    return adj, g
end



function KRegChainMatrix(net::Network)
    # connections::Array{Float64, 2} = zeros(Float64, net.nPopulations, net.nPopulations)
    # for i in 1:net.nPopulations
    #     for j in max(1, i - net.k_bar ÷ 2):min(net.nPopulations, i + net.k_bar ÷ 2)
    #         if i != j
    #             connections[i, j] = 1
    #         end
    #     end
    # end
    g = path_graph(net.nPopulations)
    adj = convert(Matrix{Float64}, adjacency_matrix(g))

    return adj, g
end

function directedPath(net::Network)
    g = path_digraph(net.nPopulations)
    adj = convert(Matrix{Float64}, adjacency_matrix(g))

    return adj, g
end


function KRegRingMatrix(net::Network)
    connections::Array{Float64, 2} = zeros(Float64, net.nPopulations, net.nPopulations)
    for i in 1:net.nPopulations
        columns = [mod(i+j-1, net.nPopulations)+1 for j in -net.k_bar÷2:net.k_bar÷2 if mod(i+j-1, net.nPopulations)+1 != i]
        connections[i, columns] .= 1
    end
    return connections
end

function updateNetwork!(populations, net, meta,s)
    epi=s.epi
    sim=s.sim
    g= net.graph

    popsRoC = Array{PopulationRoC, 1}(undef, length(populations))
    for _ in 1:sim.nTimeSteps
        # populations[1].I=.00003
        totalInfectedFlow = 0
        for (popInd, pop) in enumerate(populations)
            for (connPopInd, connWeight) in enumerate(net.connections[popInd,:])
                connPop = populations[connPopInd]
                meta.infectedFlows[popInd,connPopInd] =  connWeight * epi.μ *  (1-connPop.restrictions[pop.index]) * (1-connPop.restrictions[connPopInd])
                totalInfectedFlow += meta.infectedFlows[popInd,connPopInd]
            end
        end
        for (i, population) in enumerate(populations)
            # updatePopulation!(populations[i], connections[i,:], populations_copy, epi)
            popsRoC[i] = getPopulationRoC(population, net.connections[:,i], populations, epi,meta)

            population.S+=popsRoC[i].dS*1/sim.nTimeSteps
            population.I+=popsRoC[i].dI*1/sim.nTimeSteps
            population.R+=popsRoC[i].dR*1/sim.nTimeSteps
            population.restrictions+=popsRoC[i].restrictionsRoC*1/sim.nTimeSteps
            clamp!(populations[i].restrictions,0.0,1.0) # NEVER change these values
        end
    end
end

function computePOConnectivityHistory(restrictionsHistory,infectedHistory,populations,connections,POrder)
    # POrder = 0, 1, or 2
    nTimeSteps, nPopulations, nMobility = size(restrictionsHistory)
    POConnectivity = zeros(nTimeSteps, nPopulations)
    AveragePOConnectivity = zeros(nTimeSteps)
    for t in 1:nTimeSteps
        for i in 1:nPopulations
            for j in 1:nMobility
                POConnectivity[t, i] += restrictionsHistory[t, i, j]*restrictionsHistory[t, j, i]*connections[i,j]
                if POrder > 0 POConnectivity[t, i] *= (populations[j].size - infectedHistory[t,i]) end 
                if POrder > 1 POConnectivity[t, i] *= (populations[i].size-infectedHistory[t,j]) end
            end
            AveragePOConnectivity[t]+=POConnectivity[t, i]
        end
    end
    AveragePOConnectivity = AveragePOConnectivity./nPopulations
    return POConnectivity, AveragePOConnectivity
end
