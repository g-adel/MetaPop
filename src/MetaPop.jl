using Revise, Plots
module MetaPop
using Graphs, SparseArrays, Luxor, Karnak, Images, Plots, Plots.PlotMeasures, Printf, PrettyTables
using LinearAlgebra

# include("epidemics.jl");    include("populations.jl");  include("network.jl");
# include("visualization.jl");include("strategies.jl");   include("simulation.jl");
# include("analysis.jl");     include("multiAnalysis.jl");include("util.jl");
# include("plots.jl");        include("plotsEnsemble.jl");include("output.jl");       
# include("tests.jl");
files = ["analysis.jl", "multiAnalysis.jl", "plots.jl", "plotsEnsemble.jl", "output.jl",
            "tests.jl", "epidemics.jl", "populations.jl", "network.jl", "visualization.jl","strategies.jl", "simulation.jl"]

for file in files
    include(file)
end


mutable struct Scenario # The complete configuration of constant paramters of both system and numerical methods.
    epi::SIRS_epidemic
    net::Network
    strat::Strat
    sim::Sim
end


"""
    defineMeta()

Initialize a Metapopulation model with SIRS epidemic dynamics with specified parameters.

Creates a metapopulation system with the following components:
- SIRS epidemic model with specified transmission (β), recovery (γ), reinfection (σ) and mobility (μ) rates
- Path graph network with 10 populations and average degree 3
- Individual proportional restriction mobility strategy
- Simulation parameters including timestep and initial conditions

Returns:
- meta: Initialized Metapopulation object
- S: Scenario object containing model parameters
"""

function defineMeta()
    epi = SIRS_epidemic(β = 0.25,γ = 0.02, σ = .0, μ = 0.01)
    net = Network(nPopulations = 10, k_bar = 3, topology = PathGraph)
    strat = Strat(λ =0, mobBias = 2,strategy = IndivPropRestriction)
    sim = Sim(h = 0.1,min_h=1e-2, nDays = 1000, I₀=1e-5, critRange = 0)
    S = Scenario(epi, net, strat, sim)
    meta = Metapopulation(S = S, populations=Array{Population, 1}(undef, net.nPopulations),
                            mobilityRates = Graphs.LinAlg.adjacency_matrix(net.graph,Float64)*epi.μ, day = 1)

    return meta, S
end

"""
Instantiation, simulation, and analysis of a single epidemic Scenario
"""

function singleCaseMain()
    meta, S = defineMeta()
    initializePopulations!(meta)
    @time metaHist = simulateSystem(meta)
    data = dataAnalytics(metaHist,S)
    generatePrettyTable(data)

    combinedPlot = plotCase(data,S,meta;save=true,adaptive=S.strat.λ>0)
    return combinedPlot, data, meta
end

"""
Instantiation, simulation, and analysis of an ensemble of scenarios with arrays of values for some parameters
"""
function multiCaseMain()
    meta, S = defineMeta()
    Ss = multiScenario_λ(S)
    datas = multiSimulation1D(Ss)
    # Ss = multiScenario_μβ(S)
    # datas = multiSimulation2D(Ss)

    combinedPlot = plotMulti(datas, Ss,meta;save=false)
    return combinedPlot, datas
end

"""
Instantiation, simulation, and analysis of a single epidemic scenario with randomized network structures
"""
function ensembleCaseMain()
    meta, S = defineMeta()
    # Ss = multiScenario_strategies(S)
    Ss = multiScenario_1D(S)
    datas = multiSimulation1D(Ss)
    # Ss = multiScenario_μβ(S)
    # datas = multiSimulation2D(Ss)

    combinedPlot = plotEnsemble(datas, Ss,meta;save=true)
    return combinedPlot, datas
end

export multiCaseMain, singleCaseMain, ensembleCaseMain

end

println("¡Hola!")

@time combinedPlot, data, meta = MetaPop.singleCaseMain()

# @time combinedPlot, datas = MetaPop.multiCaseMain()

# @time combinedPlot, datas = MetaPop.ensembleCaseMain()

plot(combinedPlot)