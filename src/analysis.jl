function dataAnalytics(metaHist,S)
    nPts = length(metaHist)
    nPopulations= S.net.nPopulations
    days=zeros(Float64,nPts)
    susceptibleHistory = zeros(Float64,nPts, nPopulations)
    infectedHistory = zeros(Float64,nPts, nPopulations)
    recoveredHistory = zeros(Float64,nPts, nPopulations)
    ρsHistory = zeros(Float64,nPts, nPopulations, nPopulations)
    cumuInfMob = spzeros(Float64,nPopulations,nPopulations)
    t=1
    for meta in metaHist
        for i in 1:nPopulations
            days[t] = meta.day
            susceptibleHistory[t, i] =  meta.populations[i].S
            infectedHistory[t, i] = meta.populations[i].I
            recoveredHistory[t, i] = meta.populations[i].R
            ρsHistory[t, i, :] = meta.populations[i].ρs
            cumuInfMob[i,:] .+= meta.mobilityRates[i,:].*(meta.populations[i].I .- [p.I for p in meta.populations])
        end
        t+=1
    end
    # plot!(p, xAxis, ρsHistory[:,i+1,i], label="P($i,$(i+1))", color=colors[color_index], ylim=(0,1))
    data = Dict()
    data["metaHist"]=metaHist;    data["days"] = days; 
    data["susceptibleHistory"] = susceptibleHistory[1:nPts,:];    data["infectedHistory"] = infectedHistory[1:nPts,:];
    data["recoveredHistory"] = recoveredHistory[1:nPts,:];    data["ρsHistory"] = ρsHistory[1:nPts,:,:]; 
    data["cumuInfMob"] = cumuInfMob;
    data["downstream_ρs"] = hcat([ρsHistory[:,i+1,i] for i in 1:nPopulations-1]...)
    downstream_flows = zeros(size(data["downstream_ρs"]))
    for i in 1:nPopulations-1
        downstream_flows[:,i] = (1 .- data["downstream_ρs"][:,i]) .* infectedHistory[:,i]
    end
    data["downstream_flows"] = downstream_flows
    g = S.net.graph
    nTimeSteps, nPopulations = size(data["infectedHistory"])
    data["infectedAvgHistory"] = sum(data["infectedHistory"], dims=2)./nPopulations # TODO make it account for het. sizes
    data["susceptibleAvgHistory"] = sum(data["susceptibleHistory"], dims=2)./nPopulations
    data["recoveredAvgHistory"] = sum(data["recoveredHistory"], dims=2)./nPopulations
    # data["ρsAvgHistory"] = dropdims(sum(data["ρsHistory"],dims=3),dims=3)./sum(S.net.connections,dims=2)
    ρsAvg = zeros(nTimeSteps, nPopulations)

    ρPath = zeros(nTimeSteps, nPopulations)
    ρPath[:,1] = data["ρsHistory"][:,1,2]
    for i in 1:nPopulations-1
        ρPath[:,i+1] = data["ρsHistory"][:,i+1,i]
    end

    nodeStrength = sum(S.net.connections, dims=2)
    for t in 1:nTimeSteps
        for p in 1:nPopulations
            ρsAvg[t, p] = sum(data["ρsHistory"][t, p, :])./nodeStrength[p]
        end
    end

    data["ρsAvgHistory"] = ρsAvg
    maxInfecteds = findmax(data["infectedHistory"], dims =1)
    data["peakGlobalInfection"] = findmax(data["infectedAvgHistory"])
    data["peakInfInd"] = [index[1] for index in maxInfecteds[2]]' # reasons

    initInf = data["infectedHistory"][1, 1]

    spreadInfInd,spreadRates = find_roots(data["infectedHistory"],initInf)
    data["isolationInd"],_ = find_roots(ρPath,0.99)
    # spreadInfInd = [find_spread_index(data["infectedHistory"][:, i], initInf) for i in 1:nPopulations]
    data["spreadInfInd"] = spreadInfInd; data["1/Δt"] = spreadRates

    distances = dijkstra_shortest_paths(g, 1).dists
    data["pathLengths"] = distances
    # data["avgSpreadRate"] = 1/linear_regression_slope(distances,spreadInfInd) #TODO needs to be redefined for other networks
    data["avgSpreadRate"] = sum(spreadRates)/length(spreadRates)
    spreadRatesTotal=0
    nonNanSpreadRates = 0
    for spreadRate in spreadRates
        if !isnan(spreadRate)
            spreadRatesTotal+=spreadRate
            nonNanSpreadRates+=1
        end
    end
    data["avgSpreadRate"] = spreadRatesTotal/nonNanSpreadRates
    data["asympSpreadRate"] = 1/(spreadInfInd[end-1]-spreadInfInd[end-2])
    data["initSpreadRate"] = 1/(spreadInfInd[2]-spreadInfInd[1])
    # data["avgSpreadRate"] = 1/(spreadInfInd[end]-spreadInfInd[end-1])
    data["firstOrderSol"] = [firstOrderSol(i-1,S.epi) for i in 1:nPopulations]
    data["firstOrderApprox"] = firstOrderApprox(S.epi).*(0:nPopulations-1)
    # data["firstOrderSolSpreadRate"] = 1/linear_regression_slope(distances,data["firstOrderSol"])
    data["firstOrderSolSpreadRate"] = 1/(data["firstOrderSol"][end-1]-data["firstOrderSol"][end-2])
    data["referenceSpreadRate"] = 2.2√(S.epi.μ*(S.epi.β-S.epi.γ))
    data["second_day_restrictions"] = [ data["ρsHistory"][2,1,2] ; [data["ρsHistory"][2,i+1,i] for i=1:nPopulations-1] ]
    return data
end

using Roots, SpecialFunctions, LambertW

function firstOrderSol(j, epi)
    μ = epi.μ
    κ = epi.β - epi.γ - 2*epi.μ

    (j==0 || κ<=0) && return 0

    # Define the transformed equation to be solved for t
    # log of Sequation(t) = μ^j * t^j * exp(κ * t) / factorial(j) - 1
    equation(t) = j * log(μ) + j * log(t) + κ * t - logfactorial(j) 
    
    # Solve the equation numerically
    t_solution = find_zero(equation, 1.0)  # Initial guess is 1.0
    
    is_correct = abs(equation(t_solution)) < 1e-3
    !is_correct && error("The solution is incorrect")
    
    return t_solution
end

function firstOrderApprox(epi)
    μ = epi.μ
    κ = epi.β - epi.γ - 2*epi.μ
    κ <0 && return 0
    return lambertw(κ/(ℯ*μ))/κ
end

function linear_regression_slope(x, y)::Float64
    N = length(x)
    sum_x = sum(x)
    
    slope = (N * sum(x .* y) - sum_x * sum(y)) / (N * sum(x .^ 2) - sum_x^2)
    return slope
end

function find_roots(history, val)
    firstInd = [findfirst(x -> x >= val, history[:, i]) for i in 1:size(history,2)]
    # spreadInfInd = [x === nothing ? Inf64 : x for x in spreadInfInd]

    interpolatedIndices = []
    for i in 1:length(firstInd)
        if isnothing(firstInd[i])
            push!(interpolatedIndices, NaN)
            continue
        elseif firstInd[i]==1
            push!(interpolatedIndices, 1.0)
            continue
        end
        y1 = history[firstInd[i] - 1, i]
        y2 = history[firstInd[i], i]
        x1 = firstInd[i] - 1
        if y1 != y2
            x_zero = x1 + (val - y1)  / (y2 - y1)
            push!(interpolatedIndices, x_zero)
        else
            push!(interpolatedIndices, x1)
        end
    end

    rates = zeros(size(interpolatedIndices))
    rates[1]=NaN
    for i in 2:length(interpolatedIndices)
        rates[i]=1/(interpolatedIndices[i]-interpolatedIndices[i-1])
    end
    return interpolatedIndices, rates
end

