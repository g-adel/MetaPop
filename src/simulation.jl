Base.@kwdef mutable struct Sim
    nTimeSteps::Float64
    nDays::Int64
    I₀::Float64
    critRange::Int64
end

function simulateSystem(meta)
    sim= meta.S.sim
    metaHist = [meta]
    metaCritHist = [meta]
    newMeta=deepcopy(meta)

    for t in 1:sim.nDays
        newMeta = updateMetapop(newMeta,metaCritHist)
        push!(metaHist, deepcopy(newMeta))

        (sim.critRange>0 && t > sim.critRange) && break

        # Check if system reached steady state
        if t>1 && maximum([pop.S for pop in newMeta.populations])<.5
            maxRestrictionChange = 0.0
            for i in 1:length(metaHist[end].populations)
                maxRestrictionChange = max(maxRestrictionChange, maximum( 
                    abs.(metaHist[end].populations[i].ρs - metaHist[end-1].populations[i].ρs)))
            end
            maxInfectedChange = maximum(abs.([pop.I for pop in metaHist[end].populations]
                                         .- [pop.I for pop in metaHist[end-1].populations]))
            
            if max(maxRestrictionChange,maxInfectedChange)<0.001    
                break
            end
        end
    end
    @show meta.S.sim.nDays
    # metaHist; metaCritHist
    hist = sim.critRange > 0 ? metaCritHist : metaHist
    @show length(hist)
    return hist
end


function updateMetapop(meta,metaCritHist)
    newMeta=deepcopy(meta)
    sim, net = newMeta.S.sim, newMeta.S.net
    g= net.graph
    h = 1/sim.nTimeSteps

    for t in 1:sim.nTimeSteps
        k1 = metaRoC(newMeta)
        k2 = metaRoC(incrementMeta(newMeta,k1,h/2))
        k3 = metaRoC(incrementMeta(newMeta,k2,h/2))
        k4 = metaRoC(incrementMeta(newMeta,k3,h))
        η = 1/6 .*(k1 .+ 2.0*k2 .+ 2.0*k3 .+ k4)
        newMeta = incrementMeta(newMeta,η,h)
        newMeta.day+=1/sim.nTimeSteps
        newMeta.day<= meta.S.sim.critRange && push!(metaCritHist,newMeta)
    end
    return newMeta
end

function metaRoC(meta)
    populations = meta.populations; epi=meta.S.epi
    popsRoC = Array{PopulationRoC, 1}(undef, length(populations))

    for (popInd, pop) in enumerate(populations)
        for (connPopInd, connWeight) in enumerate(meta.S.net.connections[popInd,:])
            connPop = populations[connPopInd]
            meta.mobilityRates[popInd,connPopInd] =  connWeight * epi.μ *  (1-connPop.ρs[pop.index]) * (1-pop.ρs[connPopInd])
            # globalInfectedFlow += meta.mobilityRates[popInd,connPopInd] # WRONG
        end
    end
    # integration
    for (i, population) in enumerate(populations)
        globalInfectedFlow = 0
        # updatePopulation!(populations[i], connections[i,:], populations_copy, epi)
        # meta.mobilityRates need to be computed for popsRoC
        popsRoC[i] = getPopulationRoC(population, meta)
        # meta.day==1 && t<5 && i==2 && print(popsRoC[2].ρsRoC[1]," ")
    end
    return popsRoC
end

function incrementMeta(meta,popsRoC,timestep)
    newMeta = deepcopy(meta)
    for (i, population) in enumerate(newMeta.populations)
        globalInfectedFlow = 0
        # updatePopulation!(populations[i], connections[i,:], populations_copy, epi)
        # meta.mobilityRates need to be computed for popsRoC
        
        population.S = population.S + popsRoC[i].dS*timestep
        population.I = population.I + popsRoC[i].dI*timestep
        population.R = population.R + popsRoC[i].dR*timestep
        population.ρs = population.ρs .+ popsRoC[i].ρsRoC*timestep
        clamp!(newMeta.populations[i].ρs,0.0,1.0) # NEVER change these values
    end
    return newMeta
end


function metaSimulation(Ss)
    datas=Array{Dict,2}(undef,size(Ss))
    
    for i in 1:size(Ss,1)
        for j in 1:size(Ss,1)
            net=Ss[i,j].net
            newMeta= Metapopulation(S = Ss[i,j], populations=Array{Population, 1}(undef, net.nPopulations),
                         mobilityRates = spzeros(Float64, net.nPopulations, net.nPopulations))
            net.connections, net.graph = pathGraph(net;directed=false)
            initializePopulations!(newMeta)
            datas[i,j]=simulateSystem(newMeta)
            dataAnalytics!(datas[i,j],Ss[i,j])
        end
    end
    return datas
end
