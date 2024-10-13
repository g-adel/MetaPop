using Plots

function plotPlots(data, S;img_filename="")
    nPopulations = S.net.nPopulations
    plots = []
    push!(plots, plotTimeEvolution(data["infectedHistory"],data["infectedAvgHistory"]))
    push!(plots, plotRestrictions(data["ρsHistory"]))
    # push!(plots, plotAvgRestrictions(data["ρsAvgHistory"]))    
    push!(plots, plotInfectionIndices(data))
    push!(plots, plot_consecutive_infected(data["infectedHistory"],log_scale=false))
    push!(plots, plot_infect_ρ(data["infectedHistory"],data["ρsHistory"]))
    # push!(plots, plotRestrictionsGrid(data["ρsHistory"]))
    
    # Load the image and add it to the plots
    if img_filename != ""
        img = load(img_filename)
        img_plot = plot(img, seriestype=:image)
        push!(plots, (img_plot,800))
    end
    
    # Calculate the total height
    height = sum(p[2] for p in plots)
    
    # Extract the heights for the layout
    layout_heights = [p[2] for p in plots]./height
    layout_heights[end]=1-sum(layout_heights[1:end-1])
    
    l = Plots.grid(length(plots), 1, heights=layout_heights)
    combined_plot = plot([p[1] for p in plots]..., layout=l, size=(1000, height))
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
    return p, 400
end


function plotRestrictions(ρsHistory)
    nTimeSteps, nPopulations, _ = size(ρsHistory)
    p = plot(1:nTimeSteps, zeros(nTimeSteps), label="P1", legend=false)
    
    colors = palette(:jet, nPopulations * nPopulations)
    
    color_index = 1
    for i in 1:nPopulations
        for j in 1:nPopulations
            plot!(p, 1:nTimeSteps, ρsHistory[:, i, j], label="P($i,$j)", color=colors[color_index], ylim=(0,1))
            color_index += 1
        end
    end
    
    xlabel!(p, "Time (days)")
    ylabel!(p, "Mobility Restrictions")
    title!(p, "Evolution of Mobility Restrictions")
    return p, 400
end

function plotAvgRestrictions(ρsAvgHistory)
    nTimeSteps, nPopulations = size(ρsAvgHistory)
    p = plot(1:nTimeSteps, zeros(nTimeSteps), label="P1", legend=false)
    
    colors = palette(:jet, nPopulations)
    
    color_index = 1
    for i in 1:nPopulations
        plot!(p, 1:nTimeSteps, ρsAvgHistory[:, i], label="P($i)", color=colors[color_index], ylim=(0,1))
        color_index += 1
    end
    
    xlabel!(p, "Time (days)")
    ylabel!(p, "Mobility Restrictions")
    title!(p, "Evolution of Mobility Restrictions")
    return p, 600
end

function plotRestrictionsGrid(ρsHistory)
    nTimeSteps, nPopulations, _ = size(ρsHistory)
    layout = @layout [Plots.grid(nPopulations, nPopulations)]
    p = plot(layout=layout, size=(800, 800), framestyle=:none, ticks=nothing,legend=false)
    
    for i in 1:nPopulations
        for j in 1:nPopulations
            plot!(p, 1:nTimeSteps, ρsHistory[:, i, j], label="", subplot=(i-1)*nPopulations + j,ylim=(0,1))
        end
    end
    
    return p, 800
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
    return p, 800
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
    
    return p, 800
end

function plot_consecutive_infected(infectedHistory::Array{Float64,2}; log_scale::Bool=false)
    nTimeSteps, nPopulations = size(infectedHistory)
    p = plot()
    colors = palette(:jet, nPopulations)
    color_index = 1
    for i in 1:nPopulations-1
        plot!(p, infectedHistory[50:end,i], infectedHistory[50:end,i+1], label = "($i,$(i+1))", color = colors[color_index])
        color_index += 1
    end
    xlabel!(p, "Population i infected prevalence")
    ylabel!(p, "Population i+1 infected prevalence")
    title!(p, "Consecutive Infecteds prevalence")
    
    if log_scale
        plot!(p,xscale=:ln,yscale=:ln)
    end
    
    return p, 800
end

function plot_infect_ρ(infectedHistory::Array{Float64,2},ρsHistory::Array{Float64,3})
    nTimeSteps, nPopulations = size(infectedHistory)
    p=plot()
    colors = palette(:jet, nPopulations)
    color_index=1
    for i in 1:nPopulations-1
        plot!(p, infectedHistory[10:end,i],ρsHistory[10:end,i+1,i], label = "($i,$(i+1))", color=colors[color_index],xscale=:log10)
        color_index += 1
    end
    xlabel!(p, "Population i infected prevalence")
    ylabel!(p, "Population i+1 infected prevalence")
    title!(p, "Evolution of infected flow and ρ")
    return p, 800
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


