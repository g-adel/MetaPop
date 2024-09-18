using Revise
module MetaPop
using Graphs, Luxor, Images, Plots, Plots.PlotMeasures

include("epidemics.jl");include("populations.jl");include("network.jl");
include("visualization.jl");include("plots.jl");include("strategies.jl");
include("simulation.jl"); include("analysis.jl");
include("util.jl")

mutable struct Scenario
    epi::SIRS_epidemic
    net::Network
    strat::Strat
    sim::Sim
end


function main()
    print("¡Hola!")
    # epi = SIR_epidemic(γ = 0.05, β = 0.2, μ = .01)
    epi = SIRS_epidemic(γ = 0.05, β = 0.2, σ = .0, μ = .02)
    net = Network(; nPopulations = 10, k_bar = 2, connections = Array{Float64,2}(undef, 0, 0))
    strat = Strat(; λ = 30000, mobBias = 0) 
    sim = Sim(; nTimeSteps =100, nDays = 500)
    S = Scenario(epi, net, strat, sim)
    meta = Metapopulation(populations=Array{Population, 1}(undef, net.nPopulations), infectedFlows = zeros(Float64, net.nPopulations, net.nPopulations))
    # connections = KRegChainMatrix(net)
    net.connections = smallWorldMatrix(net)
    # println("SWM",connections,"chain",KRegChainMatrix(net))
    initializePopulations!(meta.populations,strat)
    susceptibleHistory = zeros(sim.nDays, net.nPopulations)
    infectedHistory = zeros(sim.nDays, net.nPopulations)
    recoveredHistory = zeros(sim.nDays, net.nPopulations)
    restrictionsHistory = zeros(sim.nDays, net.nPopulations, net.nPopulations)
    
    for t in 1:sim.nDays
        for i in 1:net.nPopulations
            susceptibleHistory[t, i] =  meta.populations[i].S
            infectedHistory[t, i] = meta.populations[i].I
            recoveredHistory[t, i] = meta.populations[i].R
            restrictionsHistory[t, i, :] = meta.populations[i].restrictions 
        end
        updateNetwork!(meta.populations, net, meta,S)
        
        if t>1 && maximum(susceptibleHistory[t,:])<0.03
            maxRestrictionChange = maximum(abs.(restrictionsHistory[t, :, :] .- restrictionsHistory[t-1, :, :]))
            if maxRestrictionChange<0.001    
                sim.nDays = t        
                break
            end
        end
    end
    data = Dict([("susceptibleHistory", susceptibleHistory[1:S.sim.nDays,:]), ("infectedHistory", infectedHistory[1:S.sim.nDays,:]),
     ("recoveredHistory", recoveredHistory[1:S.sim.nDays,:]), ("restrictionsHistory", restrictionsHistory[1:S.sim.nDays,:,:])])

    # data = Dict()
    # data["susceptibleHistory"] = susceptibleHistory;    data["infectedHistory"] = infectedHistory
    # data["recoveredHistory"] = recoveredHistory;    data["restrictionsHistory"] = restrictionsHistory

    P0Connectivity, AverageP0Connectivity = computePOConnectivityHistory(restrictionsHistory,infectedHistory,meta.populations,net.connections,0)
    # P1Connectivity, AverageP1Connectivity = computePOConnectivityHistory(restrictionsHistory,infectedHistory,populations,connections,1)
    P2Connectivity, AverageP2Connectivity = computePOConnectivityHistory(restrictionsHistory,infectedHistory,meta.populations,net.connections,2)
    
    # add to Dict
    data["P0Connectivity"] = P0Connectivity
    data["AverageP0Connectivity"] = AverageP0Connectivity
    # data["P1Connectivity"] = P1Connectivity
    # data["AverageP1Connectivity"] = AverageP1Connectivity
    data["P2Connectivity"] = P2Connectivity
    data["AverageP2Connectivity"] = AverageP2Connectivity

    img_file = ""
    img_file = drawNetworkPNG(meta.populations,net.connections,infectedHistory, susceptibleHistory, restrictionsHistory)
    # animate_network(populations,connections,infectedHistory, susceptibleHistory, recoveredHistory, restrictionsHistory)
    dataAnalytics!(data,net)
    # println("infection Analytics",infectionAnalytics(infectedHistory))
    # println("restrictions Analytics", restrictionsAnalytics(restrictionsHistory))
    # @show S
    plotPlots(data,S;img_filename=img_file)
    # DO NOT ADD ANYTHING HERE
end

export main

end

MetaPop.main()