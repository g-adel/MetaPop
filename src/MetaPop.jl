using Revise, Plots
module MetaPop
using Graphs, SparseArrays, Luxor, Karnak, Images, Plots, Plots.PlotMeasures, Printf, PrettyTables
using LinearAlgebra

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
    epi = SIRS_epidemic(β = 0.25,γ = 0.00, σ = .0, μ = 0.01)
    net = Network(nPopulations = 10, k_bar = 4, topology = PathGraph)
    strat = Strat(λ =4e-2, mobBias = 0.0,strategy = IndivLogRestriction)
    sim = Sim(h = 0.1,min_h=1e-2, nDays = 1000, I₀=1e-5, critRange = 0)
    S = Scenario(epi, net, strat, sim)
    meta = Metapopulation(S = S, populations=Array{Population, 1}(undef, net.nPopulations),
                            mobilityRates = Graphs.LinAlg.adjacency_matrix(net.graph,Float64)*epi.μ, day = 1)

    return meta, S
end

function singleCaseMain()
    meta, S = defineMeta()
    initializePopulations!(meta)
    @time metaHist = simulateSystem(meta)
    data = dataAnalytics(metaHist,S)
    generatePrettyTable(data)
    

    combinedPlot = plotCase(data,S,meta;save=false,adaptive=S.strat.λ>0)
    return combinedPlot, data, meta
end

function multiCaseMain()
    meta, S = defineMeta()
    Ss = multiScenario_λ(S)
    datas = multiSimulation1D(Ss)
    # Ss = multiScenario_μβ(S)
    # datas = multiSimulation2D(Ss)

    combinedPlot = plotMulti(datas, Ss,meta;save=true)
    return combinedPlot, datas
end

function ensembleCaseMain()
    meta, S = defineMeta()
    Ss = multiScenario_1D(S)
    datas = multiSimulation1D(Ss)
    # Ss = multiScenario_μβ(S)
    # datas = multiSimulation2D(Ss)

    combinedPlot = plotMulti(datas, Ss,meta;save=true)
    return combinedPlot, datas
end

export multiCaseMain, singleCaseMain

end

println("¡Hola!")

@time combinedPlot, data, meta = MetaPop.singleCaseMain()

# @time combinedPlot, datas = MetaPop.multiCaseMain()

# @time combinedPlot, datas = MetaPop.ensembleCaseMain()

plot(combinedPlot)