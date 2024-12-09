@enum GraphTopology begin
    PathGraph
    DiPathGraph
    SmallWorld
    BarabasiAlbert
end

mutable struct Network
    nPopulations::Int
    k_bar::Int
    graph
    topology::GraphTopology
    function Network(; nPopulations::Int, k_bar::Int, topology::GraphTopology)
        graph = SimpleDiGraph()
        net = new(nPopulations, k_bar, graph, topology)

        if topology == PathGraph
            net.graph = pathGraph(net; directed = false)
        elseif topology == DiPathGraph
            net.graph = pathGraph(net; directed = true)
        elseif topology == SmallWorld
            net.graph = smallWorldMatrix(net)
        elseif topology == BarabasiAlbert
            net.graph = baraAlbert(net)
        end

        return net
    end
end


function smallWorldMatrix(net::Network; β=.5)
    is_graph_connected = false
    g = SimpleGraph()
    while !is_graph_connected
        g = watts_strogatz(net.nPopulations, net.k_bar,β)
        is_graph_connected = Graphs.is_connected(g)
    end
    adj = convert(Matrix{Float64}, adjacency_matrix(g))
    return g
end

function baraAlbert(net::Network)
    g = barabasi_albert(net.nPopulations, net.k_bar)
    adj = convert(Matrix{Float64}, adjacency_matrix(g))
    return g
end



function pathGraph(net::Network; directed=false)
    if directed
        g = path_digraph(net.nPopulations)
        adj = convert(Matrix{Float64}, adjacency_matrix(g))' #FIX
    else
        g = path_graph(net.nPopulations)
        adj = convert(Matrix{Float64}, adjacency_matrix(g))
    end

    return g
end


# function KRegRingMatrix(net::Network)
#     connections::Array{Float64, 2} = zeros(Float64, net.nPopulations, net.nPopulations)
#     for i in 1:net.nPopulations
#         columns = [mod(i+j-1, net.nPopulations)+1 for j in -net.k_bar÷2:net.k_bar÷2 if mod(i+j-1, net.nPopulations)+1 != i]
#         connections[i, columns] .= 1
#     end
#     return connections
# end

# function computePOConnectivityHistory(ρsHistory,infectedHistory,populations,connections,POrder)
#     # POrder = 0, 1, or 2
#     nTimeSteps, nPopulations, nMobility = size(ρsHistory)
#     POConnectivity = zeros(nTimeSteps, nPopulations)
#     AveragePOConnectivity = zeros(nTimeSteps)
#     for t in 1:nTimeSteps
#         for i in 1:nPopulations
#             for j in 1:nMobility
#                 POConnectivity[t, i] += ρsHistory[t, i, j]*ρsHistory[t, j, i]*connections[i,j]
#                 if POrder > 0 POConnectivity[t, i] *= (populations[j].size - infectedHistory[t,i]) end 
#                 if POrder > 1 POConnectivity[t, i] *= (populations[i].size-infectedHistory[t,j]) end
#             end
#             AveragePOConnectivity[t]+=POConnectivity[t, i]
#         end
#     end
#     AveragePOConnectivity = AveragePOConnectivity./nPopulations
#     return POConnectivity, AveragePOConnectivity
# end
