function plotTimeEvolution(infected::Array{Float64,2})
    nTimeSteps, nPopulations = size(infected)

    infected_percent = infected 
    avg_infected_percent = sum(infected, dims=2)./nPopulations # calculate the average across populations

    p = plot(1:nTimeSteps, infected_percent[:, 1], label="Pop. 1")
    for i in 2:nPopulations
        plot!(p,1:nTimeSteps, infected_percent[:, i], label="Pop. $i")
    end

    # plot the average infection fraction
    plot!(p, 1:nTimeSteps, avg_infected_percent, label="Average", linestyle=:dash, linewidth=2)

    xlabel!(p,"Time (days)")
    ylabel!(p,"Infected Population Fraction")
    title!(p,"Evolution of Populations")
    return(p)
end


function plotTotalConnectivity(POConnectivity,AveragePOConnectivity,POrder,nTimeSteps,nPopulations)

    p = plot(1:nTimeSteps, POConnectivity[:, 1], label="Pop. 1")
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

