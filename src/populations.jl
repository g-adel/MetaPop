mutable struct Population
    S::Float64 # Susciptible population fraction
    I::Float64 # Infected population fraction
    R::Float64 # Recovered population fraction
    ρs#::Array{Float64, 1} # Mobility rates to other populations
    size::Float64 # Total population size
    index::Int # Index number of the population
end

@kwdef mutable struct Metapopulation
    S
    populations::Array{Population, 1}
    mobilityRates
    day::Float64
end

mutable struct PopulationRoC
    dS::Float64 # Susciptible RoC
    dI::Float64 # Infected RoC
    dR::Float64 # Recovered RoC
    ρsRoC#::Array{Float64, 1}
end
Base.zero(::Type{PopulationRoC}) = PopulationRoC(0.0, 0.0, 0.0, zeros(0))

function Base.:*(scalar::Float64, popsRoC::PopulationRoC)::PopulationRoC
    return PopulationRoC(
        popsRoC.dS * scalar,
        popsRoC.dI * scalar,
        popsRoC.dR * scalar,
        popsRoC.ρsRoC .* scalar
    )
end

function Base.:+(popsRoC1::PopulationRoC, popsRoC2::PopulationRoC)::PopulationRoC
    return PopulationRoC(
        popsRoC1.dS + popsRoC2.dS,
        popsRoC1.dI + popsRoC2.dI,
        popsRoC1.dR + popsRoC2.dR,
        popsRoC1.ρsRoC .+ popsRoC2.ρsRoC
    )
end

function initializePopulations!(meta)
    populations=meta.populations
    nPopulations = meta.S.net.nPopulations
    for i in 1:nPopulations
        populations[i] = Population(1., 0., 0., spzeros(nPopulations),1., i)
    end
    populations[1].I = meta.S.sim.I₀;
    populations[1].S = 1 - populations[1].I
end


function getPopulationRoC(pop::Population,meta::Metapopulation;globalInfFlow =0)
    net, populations, epi = meta.S.net, meta.populations, meta.S.epi;
    g=net.graph
    S = pop.S; I = pop.I; R = pop.R
    netFlowInfected = netFlowSusceptible = netFlowRecovered = 0
        # for connPopInd in inneighbors(g,pop.index)
        #     connPop = populations[connPopInd]
        #     finalMobilityRate = meta.mobilityRates[pop.index,connPopInd]
        #     netFlowSusceptible += finalMobilityRate * (connPop.S)
        #     netFlowInfected += finalMobilityRate * (connPop.I)
        #     netFlowRecovered += finalMobilityRate * (connPop.R)
        # end
        # for connPopInd in outneighbors(g,pop.index)
        #     connPop = populations[connPopInd]
        #     finalMobilityRate = meta.mobilityRates[connPopInd, pop.index]
        #     netFlowSusceptible += finalMobilityRate * (connPop.S- S)
        #     netFlowInfected += finalMobilityRate    * (connPop.I- I)
        #     netFlowRecovered += finalMobilityRate   * (connPop.R- R)
        # end
    for connPopInd in neighbors(g,pop.index)
        connPop = populations[connPopInd]
        finalMobilityRate = meta.mobilityRates[connPopInd, pop.index]

        netFlowSusceptible += finalMobilityRate * (connPop.S- S)
        netFlowInfected += finalMobilityRate    * (connPop.I- I)
        netFlowRecovered += finalMobilityRate   * (connPop.R- R)
    end
    dS = -epi.β*I*S + epi.σ*R + netFlowSusceptible 
    dI =  epi.β*I*S - epi.γ*I + netFlowInfected
    dR =  epi.γ*I   - epi.σ*R + netFlowRecovered
    ρsRoC = spzeros(meta.S.net.nPopulations)
    strategy = meta.S.strat.strategy
    if meta.S.strat.λ>0
        if strategy == UniformPropRestriction   ρsRoC = uniformPropRestriction(pop,meta)
        elseif strategy == GlobalDiffRestriction        ρsRoC = globalDiffRestriction(pop,meta,globalInfFlow)
        elseif strategy == IndivPropRestriction     ρsRoC = indivPropRestriction(pop,meta)
        elseif strategy == IndivLogRestriction      ρsRoC = indivLogRestriction(pop,meta)
        end
    end
    return  PopulationRoC(dS,dI,dR,ρsRoC) #struct
end