function plotPlots(data,nTimeSteps,net)
    plots = []
    push!(plots, plotTimeEvolution(data["infectedHistory"]))
    push!(plots, plotRestrictions(data["restrictionsHistory"]))
    push!(plots, plotRestrictionsGrid(data["restrictionsHistory"]))
    # push!(plots, plotTotalConnectivity(data["P0Connectivity"],data["AverageP0Connectivity"],0,nTimeSteps,net.nPopulations))
    # push!(plots, plotTotalConnectivity(data["P2Connectivity"],data["AverageP2Connectivity"],2,nTimeSteps,net.nPopulations))
    
    # for p in plots
    #     p = plot!(p, legend = false)
    # end
    height = 400 * length(plots)
    layout = @layout [a; b; c{0.4h}; d; e]
    plot(plots..., layout=layout, size = (600, height))
end

function plotTimeEvolution(infected::Array{Float64,2})
    nTimeSteps, nPopulations = size(infected)

    infected_percent = infected 
    avg_infected_percent = sum(infected, dims=2)./nPopulations # calculate the average across populations
    colors = palette(:jet, nPopulations)

    color_index = 1
    p = plot()
    for i in 1:nPopulations
        plot!(p,1:nTimeSteps, infected_percent[:, i], label="Pop. $i", color=colors[color_index])
        color_index += 1
    end

    # plot the average infection fraction
    plot!(p, 1:nTimeSteps, avg_infected_percent, label="Average", linestyle=:dash, linewidth=2, legendfontsize=6)

    xlabel!(p,"Time (days)")
    ylabel!(p,"Infected Population Fraction")
    title!(p,"Prevalence of Infected")
    return(p)
end

function plotRestrictions(restrictionsHistory)
    nTimeSteps, nPopulations, _ = size(restrictionsHistory)
    p = plot(1:nTimeSteps, zeros(nTimeSteps), label="P1", legend=false)
    
    colors = palette(:jet, nPopulations * nPopulations)
    
    color_index = 1
    for i in 1:nPopulations
        for j in 1:nPopulations
            plot!(p, 1:nTimeSteps, restrictionsHistory[:, i, j], label="P($i,$j)", color=colors[color_index], ylim=(0,1))
            color_index += 1
        end
    end
    
    xlabel!(p, "Time (days)")
    ylabel!(p, "Mobility Restrictions")
    title!(p, "Evolution of Mobility Restrictions")
    return p
end

function plotRestrictionsGrid(restrictionsHistory)
    nTimeSteps, nPopulations, _ = size(restrictionsHistory)
    layout = @layout [grid(nPopulations, nPopulations)]
    p = plot(layout=layout, size=(800, 800), framestyle=:none, ticks=nothing,legend=false)
    
    for i in 1:nPopulations
        for j in 1:nPopulations
            plot!(p, 1:nTimeSteps, restrictionsHistory[:, i, j], label="", subplot=(i-1)*nPopulations + j)
        end
    end
    
    return p
end

function plotTotalConnectivity(POConnectivity,AveragePOConnectivity,POrder,nTimeSteps,nPopulations)

    p = plot(1:nTimeSteps, POConnectivity[:, 1], label="Pop. 1",legend=:topright, legendbg=:transparent)
    for i in 2:nPopulations
        plot!(p, 1:nTimeSteps, POConnectivity[:, i], label="Pop. $i")
    end
    plot!(p, 1:nTimeSteps, AveragePOConnectivity.*ones(nTimeSteps), label="P"*string(POrder)*"Avg. Conn.", linestyle=:dash, linewidth=2) 
    xlabel!(p, "Time (days)")
    ylabel!(p, "Total Connectivity")
    
    (POrder == 0) && title!(p, "Evolution of M Connectivity")
    (POrder == 2) && title!(p, "Evolution of P Connectivity")
    p = plot!(p, legend=:topright)
    current_ylims = ylims(p)
    plot!(p, ylims=(0, max(current_ylims...)*1.1))
    # print ylims of plot

    return(p)
end

