Base.@kwdef mutable struct Sim
    nTimeSteps::Float64
    nDays::Int64
end

function simulateSystem(meta)
    sim, net= meta.S.sim, meta.S.net
    susceptibleHistory = zeros(sim.nDays, net.nPopulations)
    infectedHistory = zeros(sim.nDays, net.nPopulations)
    recoveredHistory = zeros(sim.nDays, net.nPopulations)
    restrictionsHistory = zeros(sim.nDays, net.nPopulations, net.nPopulations)
    
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