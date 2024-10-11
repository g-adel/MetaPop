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
    restrictionsHistory = zeros(Float64,sim.nDays, net.nPopulations, net.nPopulations)
    
    for t in 1:meta.S.sim.nDays
        for i in 1:net.nPopulations
            susceptibleHistory[t, i] =  meta.populations[i].S
            infectedHistory[t, i] = meta.populations[i].I
            recoveredHistory[t, i] = meta.populations[i].R
            restrictionsHistory[t, i, :] = meta.populations[i].restrictions 
        end
        updateNetwork!(meta)
        
        if t>1 && maximum(susceptibleHistory[t,:])<.5
            maxRestrictionChange = maximum(abs.(restrictionsHistory[t, :, :] .- restrictionsHistory[t-1, :, :]))
            maxInfectedChange = maximum(abs.(infectedHistory[t,:] .- infectedHistory[t-1,:]))
            if max(maxRestrictionChange,maxInfectedChange)<0.001    
                meta.S.sim.nDays = t        
                break
            end
        end
    end
        
    data = Dict()
    data["susceptibleHistory"] = susceptibleHistory[1:sim.nDays,:];    data["infectedHistory"] = infectedHistory[1:sim.nDays,:];
    data["recoveredHistory"] = recoveredHistory[1:sim.nDays,:];    data["restrictionsHistory"] = restrictionsHistory[1:sim.nDays,:,:]

    return data
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
