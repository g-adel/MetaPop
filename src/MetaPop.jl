using Revise, Plots
module MetaPop
using Graphs, SparseArrays, Luxor, Karnak, Images, Plots, Plots.PlotMeasures, Printf, PrettyTables

include("epidemics.jl");    include("populations.jl");  include("network.jl");
include("visualization.jl");include("strategies.jl");   include("simulation.jl");
include("analysis.jl");     include("metaAnalysis.jl"); include("util.jl");
include("plots.jl");        include("output.jl");       include("tests.jl");

mutable struct Scenario
    epi::SIRS_epidemic
    net::Network
    strat::Strat
    sim::Sim
end

function defineMeta()
    epi = SIRS_epidemic(β = 0.25,γ = 0.0, σ = .0, μ = 1/50)
    net = Network(; nPopulations = 10, k_bar = 2, topology = PathGraph)
    strat = Strat(; λ = 1e5, mobBias = 0.0,strategy = IndivDiffRestriction)
    sim = Sim(; nTimeSteps =200, nDays = 500, I₀=1e-5, critRange = 0)
    S = Scenario(epi, net, strat, sim)
    meta = Metapopulation(S = S, populations=Array{Population, 1}(undef, net.nPopulations),
                         mobilityRates = sparse(net.connections)*epi.μ, day = 1)

    return meta, S
end


function main()
    println("¡Hola!")
    meta, S = defineMeta()

    initializePopulations!(meta)
    metaHist = simulateSystem(meta)
    data = dataAnalytics(metaHist,S)
    generatePrettyTable(data)

    # Ss = multiScenario_μβ(S)
    # datas = metaSimulation(Ss)

    # plot = plot_spread_rates(datas, Ss)
    # return plot[1], datas

    img_file = ""
    # img_file = drawNetworkPNG(meta.populations,net.connections,infectedHistory, susceptibleHistory, ρsHistory)
    # img_file = drawNetworkKarnak(meta, net, data)
    # animate_network(meta.populations,meta.connections,infectedHistory, susceptibleHistory, recoveredHistory, ρsHistory)
    # println("infection Analytics",infectionAnalytics(infectedHistory))
    # println("ρs Analytics", ρsAnalytics(ρsHistory))
    # @show data
    combinedPlot = plotCase(data,S;img_filename=img_file)
    return combinedPlot, meta, data
end

export main

end
# combinedPlot, datas=MetaPop.main()
combinedPlot, meta, data = MetaPop.main()
plot(combinedPlot)