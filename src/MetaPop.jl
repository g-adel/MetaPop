using Revise, Plots
module MetaPop
using Graphs, Luxor, Karnak, Images, Plots, Plots.PlotMeasures

include("epidemics.jl");include("populations.jl");include("network.jl");
include("visualization.jl");include("plots.jl");include("strategies.jl");
include("simulation.jl"); include("analysis.jl");
include("util.jl"); include("tests.jl")

mutable struct Scenario
    epi::SIRS_epidemic
    net::Network
    strat::Strat
    sim::Sim
end


function main()
    println("¡Hola!")
    # epi = SIR_epidemic(γ = 0.05, β = 0.2, μ = .01)
    epi = SIRS_epidemic(β = 0.2,γ = 0.05, σ = .0, μ = .02)
    net = Network(; nPopulations = 10, k_bar = 2, connections = Array{Float64,2}(undef, 0, 0), graph = SimpleDiGraph())
    strat = Strat(; λ = 0000, mobBias = 0)
    sim = Sim(; nTimeSteps =50, nDays = 500)
    S = Scenario(epi, net, strat, sim) #TODO add S to meta
    meta = Metapopulation(S = S, populations=Array{Population, 1}(undef, net.nPopulations),
                         infectedFlows = zeros(Float64, net.nPopulations, net.nPopulations))
    # net.connections, net.graph = KRegChainMatrix(net)
    net.connections, net.graph = directedPath(net)
    # net.connections, net.graph = smallWorldMatrix(net)
    # net.connections, net.graph = baraAlbert(net)
    # println("SWM",connections,"chain",KRegChainMatrix(net))
    debugPlots = []
    debugVars = Dict()
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

    dataAnalytics!(data,net)

    img_file = ""
    # img_file = drawNetworkPNG(meta.populations,net.connections,infectedHistory, susceptibleHistory, restrictionsHistory)
    img_file = drawNetworkKarnak(meta, net, data)
    # animate_network(meta.populations,meta.connections,infectedHistory, susceptibleHistory, recoveredHistory, restrictionsHistory)
    # println("infection Analytics",infectionAnalytics(infectedHistory))
    # println("restrictions Analytics", restrictionsAnalytics(restrictionsHistory))
    # @show data
    combinedPlot = plotPlots(data,S;img_filename=img_file)
    return combinedPlot, meta, data
end

export main

end

combinedPlot, meta, data = MetaPop.main()

plot(combinedPlot)