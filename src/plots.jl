using Plots

function plotPlots(data, S;img_filename="")
    nPopulations = S.net.nPopulations
    plots = []
    push!(plots, plotTimeEvolution(data["infectedHistory"],data["infectedAvgHistory"]))
    push!(plots, plotRestrictions(data["restrictionsHistory"]))
    # push!(plots, plotAvgRestrictions(data["restrictionsAvgHistory"]))    
    push!(plots, plotInfectionIndices(data))
    # push!(plots, plotRestrictionsGrid(data["restrictionsHistory"]))
    
    # Load the image and add it to the plots
    if img_filename != ""
        img = load(img_filename)
        img_plot = plot(img, seriestype=:image)
        push!(plots, img_plot)
    end
    # Adjust the layout to accommodate the new plot
    height = 800 * length(plots)
    layout = @layout [a; b; c{0.4h}; d{0.4h};]
    
    combined_plot = plot(plots..., layout=layout, size=(1000, height))
    return combined_plot
end

function plotTimeEvolution(timeseries::Array{Float64,2},avg_infected_percent)
    nTimeSteps, nPopulations = size(timeseries)

    infected_percent = timeseries 
    colors = palette(:jet, nPopulations)

    color_index = 1
    p = plot()
    for i in 1:nPopulations
        plot!(p,1:nTimeSteps, infected_percent[:, i], label="Pop. $i", color=colors[color_index])
        color_index += 1
    end
    plot!(p, 1:nTimeSteps, avg_infected_percent, label="Average",
     linestyle=:dash, linewidth=2, legendfontsize=6, legend=false)
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

function plotAvgRestrictions(restrictionsAvgHistory)
    nTimeSteps, nPopulations = size(restrictionsAvgHistory)
    p = plot(1:nTimeSteps, zeros(nTimeSteps), label="P1", legend=false)
    
    colors = palette(:jet, nPopulations)
    
    color_index = 1
    for i in 1:nPopulations
        plot!(p, 1:nTimeSteps, restrictionsAvgHistory[:, i], label="P($i)", color=colors[color_index], ylim=(0,1))
        color_index += 1
    end
    
    xlabel!(p, "Time (days)")
    ylabel!(p, "Mobility Restrictions")
    title!(p, "Evolution of Mobility Restrictions")
    return p
end

function plotRestrictionsGrid(restrictionsHistory)
    nTimeSteps, nPopulations, _ = size(restrictionsHistory)
    layout = @layout [Plots.grid(nPopulations, nPopulations)]
    p = plot(layout=layout, size=(800, 800), framestyle=:none, ticks=nothing,legend=false)
    
    for i in 1:nPopulations
        for j in 1:nPopulations
            plot!(p, 1:nTimeSteps, restrictionsHistory[:, i, j], label="", subplot=(i-1)*nPopulations + j,ylim=(0,1))
        end
    end
    
    return p
end

function plotInfectionIndices(data)
    pathLengths = data["pathLengths"]
    spreadInfInd = data["spreadInfInd"]
    peakInfInd = data["peakInfInd"]

    # Create a scatter plot
    p = scatter(pathLengths, spreadInfInd, label="Spread Infection Day",
                xlabel="Path Length", ylabel="Indices", size=(1000,800),
                legendfontsize=15, tickfontsize=15, guidefontsize=15, left_margin=10Plots.mm,
                mc=:blue)
    
    scatter!(p, pathLengths, peakInfInd, label="Peak Infection Day", mc=:red)
    scatter!(p, pathLengths, data["firstOrderApprox"], label="First Order Approx.", mc=:yellow)
    scatter!(p, pathLengths[2:end], spreadInfInd[2:end].-spreadInfInd[1:end-1], label="Infection rate")
    p=plot(p,legend=false)
    return p
end

function plot_spread_rates(data::Array{Dict, 2},Ss)
    nRows, nCols = size(data)
    
    # Initialize arrays to store the values
    referenceSpreadRates = [[] for _ in 1:nCols]
    firstOrderSpreadRates = [[] for _ in 1:nCols]
    spreadRates = [[] for _ in 1:nCols]
    βs=[]
    μs=[]
    
    # Extract values from each dictionary in the 2D array
    for i in 1:nRows
        for j in 1:nCols
            push!(referenceSpreadRates[j], data[i, j]["referenceSpreadRate"])
            push!(firstOrderSpreadRates[j], data[i, j]["firstOrderSpreadRate"])
            push!(spreadRates[j], data[i, j]["spreadRate"])
            i==1 && push!(βs,Ss[1,j].epi.β)
        end
        push!(μs,Ss[i,1].epi.μ) 
    end
    @show βs μs
    # Create the plot
    p = plot(xlabel="μ", ylabel="Rate (nodes/day)", title="Spread Rates")
    
    for j in 1:nCols
        plot!(p, μs, spreadRates[j], label="Actual Rate (β= $(βs[j]))", lw=2, lc=:green, legendfontsize=6)
        scatter!(p, μs, referenceSpreadRates[j], label="Ref Rate (β= $(βs[j]))", lw=2, mc=:black)
        scatter!(p, μs, firstOrderSpreadRates[j], label="My Rate (β= $(βs[j]))", lw=2, mc=:white)
    end
    
    return p
end

# function plotTotalConnectivity(POConnectivity,AveragePOConnectivity,POrder,nTimeSteps,nPopulations)

#     p = plot(1:nTimeSteps, POConnectivity[:, 1], label="Pop. 1",legend=:topright, legendbg=:transparent)
#     for i in 2:nPopulations
#         plot!(p, 1:nTimeSteps, POConnectivity[:, i], label="Pop. $i")
#     end
#     plot!(p, 1:nTimeSteps, AveragePOConnectivity.*ones(nTimeSteps), label="P"*string(POrder)*"Avg. Conn.", linestyle=:dash, linewidth=2) 
#     xlabel!(p, "Time (days)")
#     ylabel!(p, "Total Connectivity")
    
#     (POrder == 0) && title!(p, "Evolution of M Connectivity")
#     (POrder == 2) && title!(p, "Evolution of P Connectivity")
#     p = plot!(p, legend=:topright)
#     current_ylims = ylims(p)
#     plot!(p, ylims=(0, max(current_ylims...)*1.1))
#     # print ylims of plot

#     return(p)
# end


