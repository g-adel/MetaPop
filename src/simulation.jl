Base.@kwdef mutable struct Sim
    h::Float64
    min_h::Float64
    nDays::Int64
    I₀::Float64
    critRange::Float64
end

function simulateSystem(meta)
    sim= meta.S.sim
    metaHist = [meta]
    metaCritHist = [meta]
    newMeta=deepcopy(meta)
    h = sim.h
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
    k = Array{Vector{PopulationRoC}, 1}(undef, 4)
    m = Array{Metapopulation, 1}(undef,4)
    c = [0,1/2,1/2,1]
    b = [1/6,1/3,1/3,1/6]
    a = [   0    0   0 0;
            .5   0   0 0;
            0    .5  0 0;
            0    0   1 0]
    
    t = 0
    while t < 1
        h_new = h
        # h_new=min(h,h_new/(h_scale^1))
        while true
            # println("attempt")
            m[1] = newMeta
            k[1] = metaRoC(m[1])
            m[2], OOB1 = incrementMeta(m[1],k[1],a[2,1]*h_new)
            k[2] = metaRoC(m[2])
            m[3], OOB2 = incrementMeta(m[1],k[2],a[3,2]*h_new)
            k[3] = metaRoC(m[3])
            m[4], OOB3 = incrementMeta(m[1],k[3],a[4,3]*h_new)
            k[4] = metaRoC(m[4])
            η = sum(b .* k)
            iterMeta, OOB4 = incrementMeta(newMeta,η,h_new)
            OOB = OOB1||OOB2||OOB3||OOB4
            # @show OOB1 OOB2 OOB3 OOB4 h_new
            if OOB && h_new >= meta.S.sim.min_h
                h_new*=h_scale
            else   
                newMeta=iterMeta 
                break   
            end
        end
        # iterMeta.day<2 && println("OOB: ", OOB ," ρ21: ", iterMeta.populations[2].ρs[1], ", h: ", h_new, ", t: ", t)
        # (h_new<h) && @printf("h: %.0f, d+t: %.9f, ρ21: %.9f, ρ32: %.9f \n "
        #                     , log10(h_new), newMeta.day+t,newMeta.populations[2].ρs[1],newMeta.populations[3].ρs[2])
        # h_new<h && print(h_new," ")
        # OOB && @warn("solution OOB! at day $(newMeta.day), t=$t, h=$h_new")
        newMeta.day+=h_new
        t+=h_new
        newMeta.day<= meta.S.sim.critRange && push!(metaCritHist,newMeta)
    end
    return newMeta
end

function updateMetapopRK45(meta, metaCritHist, h)
    newMeta = deepcopy(meta)
    h_new = h
    t = 0.0
    min_h = meta.S.sim.min_h
    max_h = h
    tol = 1e-6  # Error tolerance

    # Dormand-Prince coefficients
    c = [0.0, 1/5, 3/10, 4/5, 8/9, 1.0]
    a = [
        [],
        [1/5],
        [3/40, 9/40],
        [44/45, -56/15, 32/9],
        [19372/6561, -25360/2187, 64448/6561, -212/729],
        [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656]
    ]
    b = [35/384, 0.0, 500/1113, 125/192, -2187/6784, 11/84]
    b_star = [5179/57600, 0.0, 7571/16695, 393/640, -92097/339200, 187/2100]

    while t < 1.0
        h_new = min(h_new, 1.0 - t)  # Prevent overshooting
        OOB = false

        # Initialize k and m arrays
        k = Vector{PopulationRoC}[]
        m = Metapopulation[]

        # Stage 1
        m_1 = newMeta
        k_1 = metaRoC(m_1)
        push!(k, k_1)
        push!(m, m_1)

        # Stages 2 to 6
        for i in 2:6
            sum_a_k = zero(k_1)
            for j in 1:i-1
                sum_a_k += a[i][j] * k[j]
            end
            temp_meta, OOB_stage = incrementMeta(newMeta, sum_a_k, h_new)
            OOB = OOB || OOB_stage
            temp_k = metaRoC(temp_meta)
            push!(k, temp_k)
            push!(m, temp_meta)
        end

        # Compute high-order and low-order increments
        η = zero(k_1)
        η_star = zero(k_1)
        for i in 1:6
            η += b[i] * k[i]
            η_star += b_star[i] * k[i]
        end

        # Compute the next Metapopulations
        iterMeta, OOB_iter = incrementMeta(newMeta, η, h_new)
        iterMeta_star, OOB_star = incrementMeta(newMeta, η_star, h_new)
        OOB = OOB || OOB_iter || OOB_star

        # Estimate error
        error = norm(η - η_star)

        # Adaptive step size control
        if error <= tol || h_new <= min_h
            # Accept the step
            newMeta = iterMeta
            newMeta.day += h_new
            t += h_new
            if newMeta.day <= meta.S.sim.critRange
                push!(metaCritHist, newMeta)
            end

            # Update step size
            if error == 0.0
                s = 2.0
            else
                s = 0.9 * (tol / error)^(1/5)
            end
            h_new = clamp(s * h_new, min_h, max_h)
        else
            # Reject the step and reduce h_new
            h_new = max(h_new * 0.5, min_h)
        end

        if OOB
            @warn("Solution out of bounds at time t = $t")
            break
        end
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
    populations = meta.populations; epi=meta.S.epi; g=meta.S.net.graph
    popsRoC = Array{PopulationRoC, 1}(undef, length(populations))
    globalInfectedFlow = 0
    for (popInd, pop) in enumerate(populations)
        for connPopInd in neighbors(g,popInd)
            connPop = populations[connPopInd]
            # connWeight removed!! needs weighted graphs or another sparse matrix
            meta.mobilityRates[popInd,connPopInd] =  epi.μ *  (1-connPop.ρs[popInd]) * (1-pop.ρs[connPopInd]) * connPop.size
            globalInfectedFlow += abs(meta.mobilityRates[popInd,connPopInd]*pop.I) 
        end
    end
    for (i, population) in enumerate(populations)
        # meta.mobilityRates need to be computed for popsRoC
        popsRoC[i] = getPopulationRoC(population, meta;globalInfFlow = globalInfectedFlow)
    end
    return popsRoC
end


function multiSimulation1D(Ss)
    datas = Array{Dict,1}(undef,size(Ss))
    metaHists = []
    for i in 1:size(Ss,1)
        net = Ss[i].net
        nPopulations = net.nPopulations
        newMeta = Metapopulation(S=Ss[i], populations = Array{Population,1}(undef,nPopulations),
                            mobilityRates = spzeros(Float64, nPopulations, nPopulations),day=1)
        initializePopulations!(newMeta)
        push!(metaHists,simulateSystem(newMeta))
        datas[i] = dataAnalytics(metaHists[i],Ss[i])
    end
    return datas
end

function multiSimulation2D(Ss)
    datas=Array{Dict,2}(undef,size(Ss))
    metaHists = Array{Vector{Metapopulation},2}(undef,size(Ss))
    
    for i in 1:size(Ss,1)
        for j in 1:size(Ss,1)
            net=Ss[i,j].net
            newMeta= Metapopulation(S = Ss[i,j], populations=Array{Population, 1}(undef, net.nPopulations),
                         mobilityRates = spzeros(Float64, net.nPopulations, net.nPopulations),day = 1)
            net.graph = pathGraph(net;directed=false)
            initializePopulations!(newMeta)
            metaHists[i,j]=simulateSystem(newMeta)
            datas[i,j] = dataAnalytics(metaHists[i,j],Ss[i,j])
        end
    end
    return datas
end
