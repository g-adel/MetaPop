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
    epi = SIRS_epidemic(β = 0.125,γ = 0.02, σ = .0, μ = 0.01)
    net = Network(; nPopulations = 10, k_bar = 2, topology = PathGraph)
    strat = Strat(; λ = 1e10, mobBias = 0.0,strategy = IndivPropRestriction)
    sim = Sim(; h =0.1,min_h=10e-9, nDays = 500, I₀=1e-5, critRange = 0)
    S = Scenario(epi, net, strat, sim)
    meta = Metapopulation(S = S, populations=Array{Population, 1}(undef, net.nPopulations),
                            mobilityRates = Graphs.LinAlg.adjacency_matrix(net.graph,Float64)*epi.μ, day = 1)

    return meta, S
end

function metaParamsString(meta::Metapopulation)
    S = meta.S
    epi = S.epi
    net = S.net
    strat = S.strat
    sim = S.sim

    params_str = "β$(epi.β)_γ$(epi.γ)_σ$(epi.σ)_μ$(epi.μ)_nPop$(net.nPopulations)_kbar$(net.k_bar)_λ$(strat.λ)_mobBias$(strat.mobBias)_h$(sim.h)_minh$(sim.min_h)_nDays$(sim.nDays)_I0$(sim.I₀)_critRange$(sim.critRange)"
    
    return params_str
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
    # Ss = multiScenario_λ(S)
    # datas = multiSimulation1D(Ss)
    Ss = multiScenario_μβ(S)
    datas = multiSimulation2D(Ss)

    combinedPlot = plotEnsemble(datas, Ss,meta;save=false)
    return combinedPlot, datas
end

export multiCaseMain, singleCaseMain

end

println("¡Hola!")

combinedPlot, data, meta = MetaPop.singleCaseMain()

# combinedPlot, datas = MetaPop.multiCaseMain()

plot(combinedPlot)