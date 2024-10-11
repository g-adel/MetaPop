function dataAnalytics!(data,S)
    net=S.net
    g = net.graph
    nTimeSteps, nPopulations = size(data["infectedHistory"])
    data["infectedAvgHistory"] = sum(data["infectedHistory"], dims=2)./nPopulations # TODO make it account for het. sizes
    data["susceptibleAvgHistory"] = sum(data["susceptibleHistory"], dims=2)./nPopulations
    data["recoveredAvgHistory"] = sum(data["recoveredHistory"], dims=2)./nPopulations
    # data["restrictionsAvgHistory"] = dropdims(sum(data["restrictionsHistory"],dims=3),dims=3)./sum(S.net.connections,dims=2)
    restrictionsAvg = zeros(nTimeSteps, nPopulations)
    nodeStrength = sum(S.net.connections, dims=2)
    for t in 1:nTimeSteps
        for p in 1:nPopulations
            restrictionsAvg[t, p] = sum(data["restrictionsHistory"][t, p, :])./nodeStrength[p]
        end
    end
    data["restrictionsAvgHistory"] = restrictionsAvg
    maxInfecteds = findmax(data["infectedHistory"], dims =1)
    data["peakGlobalInfection"] = findmax(data["infectedAvgHistory"])
    peakInfInd = [index[1] for index in maxInfecteds[2]]
    data["peakInfInd"] = peakInfInd' # reasons
    initInf = data["infectedHistory"][1, 1]
    spreadInfInd = [findfirst(x -> x >= initInf, data["infectedHistory"][:, i]) for i in 1:nPopulations]
    spreadInfInd = [x === nothing ? Inf64 : x for x in spreadInfInd]
    data["spreadInfInd"] = spreadInfInd
    distances = dijkstra_shortest_paths(g, 1).dists
    data["pathLengths"] = distances
    data["spreadRate"] = 1/linear_regression_slope(distances,spreadInfInd) #TODO needs to be redefined for other networks
    # data["spreadRate"] = 1/(spreadInfInd[end]-spreadInfInd[end-1])
    data["firstOrderApprox"] = [firstOrderApprox(i-1,S.epi) for i in 1:nPopulations]
    data["firstOrderSpreadRate"] = 1/linear_regression_slope(distances,data["firstOrderApprox"])
    data["referenceSpreadRate"] = 2.2√(S.epi.μ*(S.epi.β-S.epi.γ))
    @show data["spreadRate"] data["firstOrderSpreadRate"] data["referenceSpreadRate"]

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

function linear_regression_slope(x, y)::Float64
    N = length(x)
    sum_x = sum(x)
    
    slope = (N * sum(x .* y) - sum_x * sum(y)) / (N * sum(x .^ 2) - sum_x^2)
    return slope
end