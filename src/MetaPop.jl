using Revise, Plots
module MetaPop
using Graphs, Luxor, Karnak, Images, Plots, Plots.PlotMeasures

include("epidemics.jl");include("populations.jl");include("network.jl");
include("visualization.jl");include("plots.jl");include("strategies.jl");
include("simulation.jl"); include("analysis.jl"); include("metaAnalysis.jl");
include("util.jl"); include("tests.jl");

mutable struct Scenario
    epi::SIRS_epidemic
    net::Network
    strat::Strat
    sim::Sim
end


function main()
    println("¡Hola!")
    epi = SIRS_epidemic(β = 0.2,γ = 0.0, σ = .0, μ = 1/50)
    net = Network(; nPopulations = 20, k_bar = 2, topology = PathGraph)
    strat = Strat(; λ = 10e10, mobBias = 0.0,strategy = IndivDiffRestriction)
    sim = Sim(; nTimeSteps =1, nDays = 500, I₀=1e-10)
    S = Scenario(epi, net, strat, sim)
    meta = Metapopulation(S = S, populations=Array{Population, 1}(undef, net.nPopulations),
                         mobilityRates = zeros(Float64, net.nPopulations, net.nPopulations))


    initializePopulations!(meta)
    data = simulateSystem(meta)

    dataAnalytics!(data,S)

    # Ss = multiScenario(S)
    # datas = metaSimulation(Ss)

    # plot = plot_spread_rates(datas, Ss)
    # return plot, datas

    img_file = ""
    # img_file = drawNetworkPNG(meta.populations,net.connections,infectedHistory, susceptibleHistory, ρsHistory)
    # img_file = drawNetworkKarnak(meta, net, data)
    # animate_network(meta.populations,meta.connections,infectedHistory, susceptibleHistory, recoveredHistory, ρsHistory)
    # println("infection Analytics",infectionAnalytics(infectedHistory))
    # println("ρs Analytics", ρsAnalytics(ρsHistory))
    # @show data
    combinedPlot = plotPlots(data,S;img_filename=img_file)
    return combinedPlot, meta, data
end

export main

end
# combinedPlot, datas=MetaPop.main()
combinedPlot, meta, data = MetaPop.main()
plot(combinedPlot)