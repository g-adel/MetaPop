using Revise, Plots
module MetaPop
using Graphs, Luxor, Karnak, Images, Plots, Plots.PlotMeasures

include("epidemics.jl");include("populations.jl");include("network.jl");
include("visualization.jl");include("plots.jl");include("strategies.jl");
include("simulation.jl"); include("analysis.jl");
include("util.jl"); include("tests.jl");
mutable struct Scenario
    epi::SIRS_epidemic
    net::Network
    strat::Strat
    sim::Sim
end



function main()
    println("¡Hola!")
    epi = SIRS_epidemic(β = 0.2,γ = 0.00, σ = .0, μ = 1/10)
    net = Network(; nPopulations = 15, k_bar = 2, connections = Array{Float64,2}(undef, 0, 0), graph = SimpleDiGraph())
    strat = Strat(; λ = 3000000, mobBias = 0.0)
    sim = Sim(; nTimeSteps =50, nDays = 500)
    S = Scenario(epi, net, strat, sim)
    meta = Metapopulation(S = S, populations=Array{Population, 1}(undef, net.nPopulations),
                         mobilityRates = zeros(Float64, net.nPopulations, net.nPopulations))
    net.connections, net.graph = pathGraph(net;directed=false)
    # net.connections, net.graph = smallWorldMatrix(net)
    # net.connections, net.graph = baraAlbert(net)

    debugPlots = []
    debugVars = Dict()

    initializePopulations!(meta.populations,strat)
    data = simulateSystem(meta)

    dataAnalytics!(data,S)

    img_file = ""
    # img_file = drawNetworkPNG(meta.populations,net.connections,infectedHistory, susceptibleHistory, restrictionsHistory)
    # img_file = drawNetworkKarnak(meta, net, data)
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