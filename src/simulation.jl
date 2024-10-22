Base.@kwdef mutable struct Sim
    nTimeSteps::Float64
    nDays::Int64
    I₀::Float64
    critRange::Float64
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
        if t>100
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
    h = 1/newMeta.S.sim.nTimeSteps
    t=0
    k=Array{Vector{PopulationRoC}, 1}(undef, 4)
    m=Array{Metapopulation, 1}(undef,4)
    c = [0,1/2,1/2,1]
    b = [1/6,1/3,1/3,1/6]
    a = [   0    0   0 0;
            .5   0   0 0;
            0    .5  0 0;
            0    0   1 0]
    
    # while OOB2
    #     h=h/2
    #     m[2],OOB2 = incrementMeta(newMeta,k[1],c[2]*h)
    #     @show t h meta.day
    # end
    while t < 1
        m[1]=newMeta
        OOB2=true
        k[1] = metaRoC(m[1])
        m[2],OOB2 = incrementMeta(m[1],k[1],a[2,1]*h)
        k[2] = metaRoC(m[2])
        m[3] = incrementMeta(m[1],k[2],a[3,1]*h)[1]
        k[3] = metaRoC(m[3])
        m[4] = incrementMeta(m[1],k[3],a[4,1]*h)[1]
        k[4] = metaRoC(m[4])
        η = sum(b .* k)
        newMeta, OOB = incrementMeta(newMeta,η,h)
        OOB && println("solution OOB",newMeta)
        newMeta.day+=h
        t+=h
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
        # meta.mobilityRates need to be computed for popsRoC
        popsRoC[i] = getPopulationRoC(population, meta)
    end
    return popsRoC
end

function incrementMeta(meta,popsRoC,timestep)
    newMeta = deepcopy(meta)
    OOB = false # Out Of Bounds in the Urn
    for (i, pop) in enumerate(newMeta.populations)
        globalInfectedFlow = 0
        # updatePopulation!(populations[i], connections[i,:], populations_copy, epi)
        # meta.mobilityRates need to be computed for popsRoC
        pop.S = pop.S + popsRoC[i].dS*timestep
        pop.I = pop.I + popsRoC[i].dI*timestep
        pop.R = pop.R + popsRoC[i].dR*timestep
        pop.ρs = pop.ρs .+ popsRoC[i].ρsRoC*timestep
        if maximum(pop.ρs) > 1 || max(pop.S,pop.I,pop.R)> 1 || min(pop.S,pop.I,pop.R)<0
            OOB = true
        end
        clamp!(newMeta.populations[i].ρs,0.0,1.0) # NEVER change these values
    end
    return newMeta, OOB
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
