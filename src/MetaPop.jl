using Revise, Plots
module MetaPop
using Graphs, SparseArrays, Luxor, Karnak, Images, Plots, Plots.PlotMeasures, Printf, PrettyTables

include("epidemics.jl");    include("populations.jl");  include("network.jl");
include("visualization.jl");include("strategies.jl");   include("simulation.jl");
include("analysis.jl");     include("multiAnalysis.jl");include("util.jl");
include("plots.jl");        include("plotsEnsemble.jl");include("output.jl");       
include("tests.jl");

mutable struct Scenario
    epi::SIRS_epidemic
    net::Network
    strat::Strat
    sim::Sim
end

function defineMeta()
    epi = SIRS_epidemic(β = 0.25,γ = 0.0, σ = .0, μ = 1/100)
    net = Network(; nPopulations = 8, k_bar = 2, topology = PathGraph)
    strat = Strat(; λ = 1e10, mobBias = 0.0,strategy = IndivPropRestriction)
    sim = Sim(; h =0.1,min_h=10e-10, nDays = 500, I₀=1e-5, critRange = 0)
    S = Scenario(epi, net, strat, sim)
    meta = Metapopulation(S = S, populations=Array{Population, 1}(undef, net.nPopulations),
                            mobilityRates = sparse(net.connections)*epi.μ, day = 1)

    return meta, S
end

function singleCaseMain()
    meta, S = defineMeta()
    initializePopulations!(meta)
    metaHist = simulateSystem(meta)
    data = dataAnalytics(metaHist,S)
    generatePrettyTable(data)

    combinedPlot = plotCase(data,S)
    return combinedPlot, data
end

function multiCaseMain()
    _, S = defineMeta()
    Ss = multiScenario_λ(S)
    datas = multiSimulation1D(Ss)

    combinedPlot = plotEnsemble(datas, Ss)
    return combinedPlot, datas
end

export multiCaseMain, singleCaseMain


end

println("¡Hola!")

combinedPlot, data = MetaPop.singleCaseMain()
# combinedPlot, datas = MetaPop.multiCaseMain()

plot(combinedPlot)