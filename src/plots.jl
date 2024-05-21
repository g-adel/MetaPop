function plotTimeEvolution(infected::Array{Float64,2})
    nTimeSteps, nPopulations = size(infected)

    infected_percent = infected 
    p = plot(1:nTimeSteps, infected_percent[:, 1], label="Pop. 1")
    for i in 2:nPopulations
        plot!(p,1:nTimeSteps, infected_percent[:, i], label="Pop. $i")
    end
    # change xlabel
    xlabel!(p,"Time (days)")
    ylabel!(p,"Population Fraction")
    title!(p,"Time Evolution of Populations")
    display(p)
end