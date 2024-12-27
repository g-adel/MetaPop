
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
            net.graph = smallWorldMatrix(net; β=0.5)
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
    g = barabasi_albert(net.nPopulations, net.k_bar÷2)
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


function resistance_distances(g::Graphs.SimpleGraph)
    n = nv(g)
    # Build Laplacian
    L = zeros(Float64, n, n)
    for i in 1:n
        deg_i = degree(g, i)
        L[i,i] = deg_i
        for j in neighbors(g, i)
            L[i,j] = -1
        end
    end
    # Pseudoinverse
    Lpseudo = pinv(L)
    # Compute distances
    distances = Vector{Float64}(undef, n)
    for j in 1:n
        distances[j] = Lpseudo[1,1] + Lpseudo[j,j] - 2 * Lpseudo[1,j]
    end
    return distances
end

# function KRegRingMatrix(net::Network)
#     connections::Array{Float64, 2} = zeros(Float64, net.nPopulations, net.nPopulations)
#     for i in 1:net.nPopulations
#         columns = [mod(i+j-1, net.nPopulations)+1 for j in -net.k_bar÷2:net.k_bar÷2 if mod(i+j-1, net.nPopulations)+1 != i]
#         connections[i, columns] .= 1
#     end
#     return connections
# end
