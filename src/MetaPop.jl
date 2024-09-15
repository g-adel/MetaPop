using Revise
module MetaPop
using Luxor, Images, Plots, Plots.PlotMeasures

include("epidemics.jl");include("populations.jl");include("network.jl");
include("networkViz.jl");include("plots.jl");include("strategies.jl")
include("util.jl")

function main()
    # epi = SIR_epidemic(0.05, 0.2,0.01) # γ, β, μ
    epi = SIRS_epidemic(0.05, 0.2, .01,0.01)# γ, β, σ, μ
    net = Network(7,2) # nPopulations, k_bar
    strat = Strat(1000,0.0) # λ, mobBias
    populations::Array{Population, 1} = Array{Population, 1}(undef, net.nPopulations)
    infectedFlows::Array{Float64, 2} = zeros(Float64, net.nPopulations, net.nPopulations)
    nTimeSteps::Int = 200
    connections = KRegChainMatrix(net)

    initializePopulations!(populations,strat)
    susceptibleHistory = zeros(nTimeSteps, net.nPopulations)
    infectedHistory = zeros(nTimeSteps, net.nPopulations)
    recoveredHistory = zeros(nTimeSteps, net.nPopulations)
    mobilityRestrictionsHistory = zeros(nTimeSteps, net.nPopulations, net.nPopulations)
    data = Dict([("susceptibleHistory", susceptibleHistory), ("infectedHistory", infectedHistory), ("recoveredHistory", recoveredHistory), ("mobilityRestrictionsHistory", mobilityRestrictionsHistory)])
    
    for t in 1:nTimeSteps
        updateNetwork!(populations,infectedFlows, connections,epi)
        for i in 1:net.nPopulations
            susceptibleHistory[t, i] =  populations[i].S
            infectedHistory[t, i] = populations[i].I
            recoveredHistory[t, i] = populations[i].R
            mobilityRestrictionsHistory[t, i, :] = populations[i].mobilityRestrictions
        end
    end

    P0Connectivity, AverageP0Connectivity = computePOConnectivityHistory(mobilityRestrictionsHistory,infectedHistory,populations,connections,0)
    # P1Connectivity, AverageP1Connectivity = computePOConnectivityHistory(mobilityRestrictionsHistory,infectedHistory,populations,connections,1)
    P2Connectivity, AverageP2Connectivity = computePOConnectivityHistory(mobilityRestrictionsHistory,infectedHistory,populations,connections,2)
    
    # add to Dict
    data["P0Connectivity"] = P0Connectivity
    data["AverageP0Connectivity"] = AverageP0Connectivity
    # data["P1Connectivity"] = P1Connectivity
    # data["AverageP1Connectivity"] = AverageP1Connectivity
    data["P2Connectivity"] = P2Connectivity
    data["AverageP2Connectivity"] = AverageP2Connectivity


    # drawNetworkPNG(populations,connections,infectedHistory, moInbilityRatesHistory)
    # animate_network(populations,connections,infectedHistory, susceptibleHistory, recoveredHistory, mobilityRestrictionsHistory)
    print("¡Hola!")      
    plotPlots(data,nTimeSteps,net)
end

export main

end

MetaPop.main()