using Revise
module MetaPop
using Luxor, Images, Plots, Plots.PlotMeasures

include("epidemics.jl");include("populations.jl");include("network.jl");
include("visualization.jl");include("plots.jl");include("strategies.jl")
include("simulation.jl"); include("util.jl")

function main()
    # epi = SIR_epidemic(γ = 0.05, β = 0.2, μ = .01)
    epi = SIRS_epidemic(γ = 0.05, β = 0.2, σ = .0, μ = .01)
    net = Network(; nPopulations = 7, k_bar = 2)
    strat = Strat(; λ = 30000, mobBias = 3) 
    sim = Sim(; nTimeSteps = 40, nDays = 200)
    populations::Array{Population, 1} = Array{Population, 1}(undef, net.nPopulations)
    infectedFlows::Array{Float64, 2} = zeros(Float64, net.nPopulations, net.nPopulations)
    connections = KRegChainMatrix(net)
    
    initializePopulations!(populations,strat)
    susceptibleHistory = zeros(sim.nDays, net.nPopulations)
    infectedHistory = zeros(sim.nDays, net.nPopulations)
    recoveredHistory = zeros(sim.nDays, net.nPopulations)
    restrictionsHistory = zeros(sim.nDays, net.nPopulations, net.nPopulations)
    data = Dict([("susceptibleHistory", susceptibleHistory), ("infectedHistory", infectedHistory), ("recoveredHistory", recoveredHistory), ("restrictionsHistory", restrictionsHistory)])
    
    for t in 1:sim.nDays
        for i in 1:net.nPopulations
            susceptibleHistory[t, i] =  populations[i].S
            infectedHistory[t, i] = populations[i].I
            recoveredHistory[t, i] = populations[i].R
            restrictionsHistory[t, i, :] = populations[i].restrictions 
        end

        updateNetwork!(populations,infectedFlows, connections,epi,sim)
    end

    P0Connectivity, AverageP0Connectivity = computePOConnectivityHistory(restrictionsHistory,infectedHistory,populations,connections,0)
    # P1Connectivity, AverageP1Connectivity = computePOConnectivityHistory(restrictionsHistory,infectedHistory,populations,connections,1)
    P2Connectivity, AverageP2Connectivity = computePOConnectivityHistory(restrictionsHistory,infectedHistory,populations,connections,2)
    
    # add to Dict
    data["P0Connectivity"] = P0Connectivity
    data["AverageP0Connectivity"] = AverageP0Connectivity
    # data["P1Connectivity"] = P1Connectivity
    # data["AverageP1Connectivity"] = AverageP1Connectivity
    data["P2Connectivity"] = P2Connectivity
    data["AverageP2Connectivity"] = AverageP2Connectivity


    # drawNetworkPNG(populations,connections,infectedHistory, moInbilityRatesHistory)
    # animate_network(populations,connections,infectedHistory, susceptibleHistory, recoveredHistory, restrictionsHistory)
    print("¡Hola!")      
    plotPlots(data,sim.nDays,net)
end

export main

end

MetaPop.main()