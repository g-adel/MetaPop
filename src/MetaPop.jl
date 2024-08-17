using Revise
using Luxor, Images, Plots
foreach(include, ["util.jl", "epidemics.jl", 
        "populations.jl", "networkViz.jl", "network.jl", "plots.jl", 
        "strategies.jl"])

function main()
    epi = SIR_epidemic(0.05, 0.2,0.01)
    net = Network(11,4)
    populations::Array{Population, 1} = Array{Population, 1}(undef, net.nPopulations)
    connections::Array{Float64, 2} = zeros(Float64, net.nPopulations, net.nPopulations)

    nTimeSteps::Int = 200

    initializePopulations!(populations)
    fillConnectionMatrix!(connections,net.nPopulations, net.k_bar)
    susceptibleHistory = zeros(nTimeSteps, net.nPopulations)
    infectedHistory = zeros(nTimeSteps, net.nPopulations)
    recoveredHistory = zeros(nTimeSteps, net.nPopulations)
    mobilityRatesHistory = zeros(nTimeSteps, net.nPopulations, net.nPopulations)

    for t in 1:nTimeSteps
        updateNetwork!(populations, connections,epi)
        for i in 1:net.nPopulations
            susceptibleHistory[t, i] =  populations[i].S
            infectedHistory[t, i] = populations[i].I
            recoveredHistory[t, i] = populations[i].R
            mobilityRatesHistory[t, i, :] = populations[i].mobilityRates
        end
    end

    P0Connectivity, AverageP0Connectivity = computePOConnectivityHistory(mobilityRatesHistory,infectedHistory,populations,connections,0)
    # P1Connectivity, AverageP1Connectivity = computePOConnectivityHistory(mobilityRatesHistory,infectedHistory,populations,connections,1)
    P2Connectivity, AverageP2Connectivity = computePOConnectivityHistory(mobilityRatesHistory,infectedHistory,populations,connections,2)



    # drawNetworkPNG(populations,connections,infectedHistory, moInbilityRatesHistory)
    # animate_network(populations,connections,infectedHistory, susceptibleHistory, recoveredHistory, mobilityRatesHistory)

    plots = []
    push!(plots, plotTimeEvolution(infectedHistory))
    push!(plots, plotTotalConnectivity(P0Connectivity,AverageP0Connectivity,0,nTimeSteps,net.nPopulations))
    push!(plots, plotTotalConnectivity(P2Connectivity,AverageP2Connectivity,2,nTimeSteps,net.nPopulations))

    for p in plots
        p = plot!(p, legend=false)
    end
    plot(plots..., layout = (length(plots), 1), size = (500, 1000))
end

main()