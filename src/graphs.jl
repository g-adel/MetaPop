using Graphs

g = watts_strogatz(10, 4, 0.3)
adj = convert(Matrix{Float64}, adjacency_matrix(g))
adj[2,:]

for (i,v) in enumerate(adj[2,:])
    println("i",i,"v",v)
end