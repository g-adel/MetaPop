function dataAnalytics!(data,S)
    net=S.net
    g = net.graph
    nTimeSteps, nPopulations = size(data["infectedHistory"])
    data["infectedAvgHistory"] = sum(data["infectedHistory"], dims=2)./nPopulations # TODO make it account for het. sizes
    data["susceptibleAvgHistory"] = sum(data["susceptibleHistory"], dims=2)./nPopulations
    data["recoveredAvgHistory"] = sum(data["recoveredHistory"], dims=2)./nPopulations
    maxInfecteds = findmax(data["infectedHistory"], dims =1)
    peakInfInd = [index[1] for index in maxInfecteds[2]]
    data["peakInfInd"] = peakInfInd' # reasons
    initInf = data["infectedHistory"][1, 1]
    spreadInfInd = [findfirst(x -> x >= initInf, data["infectedHistory"][:, i]) for i in 1:nPopulations]
    spreadInfInd = [x === nothing ? 0 : x for x in spreadInfInd]
    data["spreadInfInd"] = spreadInfInd
    distances = dijkstra_shortest_paths(g, 1).dists
    data["pathLengths"] = distances
end