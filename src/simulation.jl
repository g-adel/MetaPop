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
    h = 1/sim.nTimeSteps
    for t in 1:sim.nDays
        newMeta = updateMetapop(newMeta,metaCritHist,h)
        push!(metaHist, deepcopy(newMeta))

        (sim.critRange>0 && t > sim.critRange) && break

        # Check if system reached steady state
        if t>100 && maximum([pop.S for pop in newMeta.populations])<0.9
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
    # metaHist; metaCritHist
    hist = sim.critRange > 0 ? metaCritHist : metaHist
    return hist
end


function updateMetapop(meta,metaCritHist, h)
    newMeta=deepcopy(meta)
    h_scale = 0.1
    h_new = h
    t=0
    k=Array{Vector{PopulationRoC}, 1}(undef, 4)
    m=Array{Metapopulation, 1}(undef,4)
    c = [0,1/2,1/2,1]
    b = [1/6,1/3,1/3,1/6]
    a = [   0    0   0 0;
            .5   0   0 0;
            0    .5  0 0;
            0    0   1 0]
    
    
    while t < 1
        OOB=OOB2=OOB3=OOB4=true
        iterMeta=undef
        firstTime = true
        h_new = h        
        # h_new=min(h,h_new/(h_scale^1))
        while firstTime || ((OOB||OOB2||OOB3||OOB4) && h_new >= h*h_scale^20)
            # println("attempt")
            m[1] = newMeta
            k[1] = metaRoC(m[1])
            m[2], OOB2 = incrementMeta(m[1],k[1],a[2,1]*h_new)
            k[2] = metaRoC(m[2])
            m[3], OOB3 = incrementMeta(m[1],k[2],a[3,1]*h_new)
            k[3] = metaRoC(m[3])
            m[4], OOB4 = incrementMeta(m[1],k[3],a[4,1]*h_new)
            k[4] = metaRoC(m[4])
            η = sum(b .* k)
            iterMeta, OOB = incrementMeta(newMeta,η,h_new)
            h_new = OOB ? h_new*h_scale : h_new
            
            firstTime = false
        end
        # iterMeta.day<2 && println("OOB: ", OOB ," ρ21: ", iterMeta.populations[2].ρs[1], ", h: ", h_new, ", t: ", t)
        newMeta=iterMeta
        (h_new<h && newMeta.day>2) && @printf("h: %.0f, t: %.9f, ρ21: %.9f, ρ32: %.9f \n "
                            , log10(h_new), t,newMeta.populations[2].ρs[1],newMeta.populations[3].ρs[2])
        # h_new<h && print(h_new," ")
        # OOB && @warn("solution OOB! at day $(newMeta.day), t=$t, h=$h_new")
        newMeta.day+=h_new
        t+=h_new
        newMeta.day<= meta.S.sim.critRange && push!(metaCritHist,newMeta)
    end
    return newMeta
end

function incrementMeta(meta,popsRoC,timestep)
    newMeta = deepcopy(meta)
    OOB = false # Out Of Bounds Urn
    for (i, pop) in enumerate(newMeta.populations)
        globalInfectedFlow = 0
        # updatePopulation!(populations[i], connections[i,:], populations_copy, epi)
        # meta.mobilityRates need to be computed for popsRoC
        pop.S = pop.S + popsRoC[i].dS*timestep
        pop.I = pop.I + popsRoC[i].dI*timestep
        pop.R = pop.R + popsRoC[i].dR*timestep
        pop.ρs = pop.ρs .+ popsRoC[i].ρsRoC*timestep
        if maximum(pop.ρs) > 1 || minimum(pop.ρs) < 0 || any(isnan, pop.ρs) ||
            max(pop.S,pop.I,pop.R)> 1 || min(pop.S,pop.I,pop.R)<0
            OOB = true
        end
        # if !OOB && maximum(pop.ρs) >= 1
        #     println("maximum ρ: ", maximum(pop.ρs), " minimim ρ: ", minimum(pop.ρs))
        # end
        clamp!(pop.ρs,0.0,1.0) # NEVER change these values
    end
    return newMeta, OOB
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


function multiSimulation(Ss)
    datas=Array{Dict,2}(undef,size(Ss))
    metaHists = Array{Metapopulation,2}(undef,size(Ss))
    
    for i in 1:size(Ss,1)
        for j in 1:size(Ss,1)
            net=Ss[i,j].net
            newMeta= Metapopulation(S = Ss[i,j], populations=Array{Population, 1}(undef, net.nPopulations),
                         mobilityRates = spzeros(Float64, net.nPopulations, net.nPopulations))
            net.connections, net.graph = pathGraph(net;directed=false)
            initializePopulations!(newMeta)
            metaHists[i,j]=simulateSystem(newMeta)
            datas[i,j] = dataAnalytics(metaHists[i,j],Ss[i,j])
        end
    end
    return datas
end

function multiSimulation1D(Ss)
    datas = Array{Dict,1}(undef,size(Ss))
    metaHists = []
    for i in 1:size(Ss,1)
        net = Ss[i].net
        nPopulations = net.nPopulations
        newMeta = Metapopulation(S=Ss[i], populations = Array{Population,1}(undef,nPopulations),
                            mobilityRates = spzeros(Float64, nPopulations, nPopulations),day=1)
        net.connections, net.graph = pathGraph(Ss[i].net;directed=false)
        initializePopulations!(newMeta)
        push!(metaHists,simulateSystem(newMeta))
        datas[i] = dataAnalytics(metaHists[i],Ss[i])
    end
    return datas
end