Base.@kwdef mutable struct Sim
    nTimeSteps::Float64
    nDays::Int64
    I₀::Float64
end

function simulateSystem(meta)
    sim, net= meta.S.sim, meta.S.net
    susceptibleHistory = zeros(Float64,sim.nDays, net.nPopulations)
    infectedHistory = zeros(Float64,sim.nDays, net.nPopulations)
    recoveredHistory = zeros(Float64,sim.nDays, net.nPopulations)
    ρsHistory = zeros(Float64,sim.nDays, net.nPopulations, net.nPopulations)
    
    for t in 1:meta.S.sim.nDays
        for i in 1:net.nPopulations
            susceptibleHistory[t, i] =  meta.populations[i].S
            infectedHistory[t, i] = meta.populations[i].I
            recoveredHistory[t, i] = meta.populations[i].R
            ρsHistory[t, i, :] = meta.populations[i].ρs 
        end
        updateMetapop!(meta)
        
        # Check if system reached steady state
        if t>1 && maximum(susceptibleHistory[t,:])<.5
            maxRestrictionChange = maximum(abs.(ρsHistory[t, :, :] .- ρsHistory[t-1, :, :]))
            maxInfectedChange = maximum(abs.(infectedHistory[t,:] .- infectedHistory[t-1,:]))
            if max(maxRestrictionChange,maxInfectedChange)<0.001    
                meta.S.sim.nDays = t        
                break
            end
        end
    end
        
    data = Dict()
    data["susceptibleHistory"] = susceptibleHistory[1:sim.nDays,:];    data["infectedHistory"] = infectedHistory[1:sim.nDays,:];
    data["recoveredHistory"] = recoveredHistory[1:sim.nDays,:];    data["ρsHistory"] = ρsHistory[1:sim.nDays,:,:]

    return data
end


function updateMetapop!(meta)
    populations = meta.populations
    epi, sim, net = meta.S.epi, meta.S.sim, meta.S.net
    g= net.graph

    popsRoC = Array{PopulationRoC, 1}(undef, length(populations))
    for _ in 1:sim.nTimeSteps
        # integration
        for (i, population) in enumerate(populations)
            globalInfectedFlow = 0
            for (popInd, pop) in enumerate(populations)
                for (connPopInd, connWeight) in enumerate(net.connections[popInd,:])
                    connPop = populations[connPopInd]
                    meta.mobilityRates[popInd,connPopInd] =  connWeight * epi.μ *  (1-connPop.ρs[pop.index]) * (1-pop.ρs[connPopInd])
                    # globalInfectedFlow += meta.mobilityRates[popInd,connPopInd] # WRONG
                end
            end
            # updatePopulation!(populations[i], connections[i,:], populations_copy, epi)
            popsRoC[i] = getPopulationRoC(population, meta)
            
            population.S+=popsRoC[i].dS*1/sim.nTimeSteps
            population.I+=popsRoC[i].dI*1/sim.nTimeSteps
            population.R+=popsRoC[i].dR*1/sim.nTimeSteps
            population.ρs = population.ρs .+ popsRoC[i].ρsRoC*1/sim.nTimeSteps
            clamp!(populations[i].ρs,0.0,1.0) # NEVER change these values
        end
    end
end



function metaSimulation(Ss)
    datas=Array{Dict,2}(undef,size(Ss))
    
    for i in 1:size(Ss,1)
        for j in 1:size(Ss,1)
            net=Ss[i,j].net
            newMeta= Metapopulation(S = Ss[i,j], populations=Array{Population, 1}(undef, net.nPopulations),
                         mobilityRates = zeros(Float64, net.nPopulations, net.nPopulations))
            net.connections, net.graph = pathGraph(net;directed=false)
            initializePopulations!(newMeta)
            datas[i,j]=simulateSystem(newMeta)
            dataAnalytics!(datas[i,j],Ss[i,j])
        end
    end
    return datas
end
