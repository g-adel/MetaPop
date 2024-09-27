function dataAnalytics!(data,S)
    net=S.net
    g = net.graph
    nTimeSteps, nPopulations = size(data["infectedHistory"])
    data["infectedAvgHistory"] = sum(data["infectedHistory"], dims=2)./nPopulations # TODO make it account for het. sizes
    data["susceptibleAvgHistory"] = sum(data["susceptibleHistory"], dims=2)./nPopulations
    data["recoveredAvgHistory"] = sum(data["recoveredHistory"], dims=2)./nPopulations
    maxInfecteds = findmax(data["infectedHistory"], dims =1)
    data["peakGlobalInfection"] = findmax(data["infectedAvgHistory"])
    peakInfInd = [index[1] for index in maxInfecteds[2]]
    data["peakInfInd"] = peakInfInd' # reasons
    initInf = data["infectedHistory"][1, 1]
    spreadInfInd = [findfirst(x -> x >= initInf, data["infectedHistory"][:, i]) for i in 1:nPopulations]
    spreadInfInd = [x === nothing ? 0 : x for x in spreadInfInd]
    data["spreadInfInd"] = spreadInfInd
    distances = dijkstra_shortest_paths(g, 1).dists
    data["pathLengths"] = distances
    data["firstOrderApprox"] = [firstOrderApprox(i-1,S.epi) for i in 1:nPopulations]
end

using Roots, SpecialFunctions

function firstOrderApprox(p, epi)
    μ = epi.μ
    κ = epi.β - epi.γ - 2*epi.μ

    (p==0 || κ<=0) && return 0

    # Define the transformed equation to be solved for t
    # log of Sequation(t) = μ^p * t^p * exp(κ * t) / factorial(p) - 1
    equation(t) = p * log(μ) + p * log(t) + κ * t - logfactorial(p) 
    
    # Solve the equation numerically
    t_solution = find_zero(equation, 1.0)  # Initial guess is 1.0
    
    # Verify the solution
    is_correct = abs(equation(t_solution)) < 1e-3
    !is_correct && error("The solution is incorrect")
    
    return t_solution
end