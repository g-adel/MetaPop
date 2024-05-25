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


function plotTotalConnectivity(mobilityRatesHistory::Array{Float64,3},connections::Array{Float64, 2})
    nTimeSteps, nPopulations, _ = size(mobilityRatesHistory)

    totalConnectivity = zeros(nTimeSteps, nPopulations)
    AverageConnectivity = zeros(nTimeSteps)
    for i in 1:nTimeSteps
        for j in 1:nPopulations
            for k in 1:size(mobilityRatesHistory, 3)
                totalConnectivity[i, j] += mobilityRatesHistory[i, j, k]*connections[j,k]
            end
            AverageConnectivity[i]+=totalConnectivity[i, j]
        end
    end
    AverageConnectivity = AverageConnectivity./nPopulations
    p = plot(1:nTimeSteps, totalConnectivity[:, 1], label="Pop. 1")
    for i in 2:nPopulations
        plot!(p, 1:nTimeSteps, totalConnectivity[:, i], label="Pop. $i")
    end
    plot!(p, 1:nTimeSteps, AverageConnectivity.*ones(nTimeSteps), label="Avg. Conn.", linestyle=:dash, linewidth=2) 
    xlabel!(p, "Time (days)")
    ylabel!(p, "Total Connectivity")
    title!(p, "Evolution of Total Connectivity")
    p = plot!(p, legend=:topright)
    
    return(p)
end