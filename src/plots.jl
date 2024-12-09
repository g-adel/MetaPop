function alignAndCombinePlots(plots)
    gr()
    height = sum(p[2] for p in plots)
    for p in plots
        plot!(p[1], left_margin=20Plots.mm, bottom_margin=10Plots.mm, tickfontsize=10, guidefontsize=15)
    end
    
    # Extract the heights for the layout
    layout_heights = [p[2] for p in plots]./height
    layout_heights[end]=1-sum(layout_heights[1:end-1])
    
    l = Plots.grid(length(plots), 1, heights=layout_heights)
    combined_plot = plot([p[1] for p in plots]..., layout=l, size=(800, height))
    return combined_plot
end

function savePlots(plots, meta)
    # Ensure the directory exists
    dir = ".\\Thesis\\Figures"
    isdir(dir) || mkpath(dir)

    # Get the meta params string
    params_str = metaParamsString(meta)

    # Save each plot with a unique name
    for (i, plot) in enumerate(plots)
        filename = joinpath(dir, "$(plot[3])_$(params_str).pdf")
        savefig(plot[1], filename)
    end
end

function plotCase(data,S,meta;save=false,adaptive=false)
    plots = []
    push!(plots, plotInfEvolution(data))
    push!(plots, plotInfEvolution(data;yLog=true))
    S.sim.critRange>0 && push!(plots, plotInfEvolution(data;yLog=true,xLog=true))
    # push!(plots, plotInfectedFlow(data))
    # push!(plots, plotAvgRestrictions(data["ρsAvgHistory"]))    
    push!(plots, plotInfectionDays(data,S.epi))
    push!(plots, plotSpreadRates(data))
    # push!(plots, plot_cumulative_flow(data))
    # push!(plots, plot_consecutive_infected(data["infectedHistory"],log_scale=false))
    if adaptive
        push!(plots, plotRestrictions(data))
        push!(plots, plot_infect_ρ(data["infectedHistory"],data["ρsHistory"]))
    end
        # push!(plots, plot_infect_flow(data["infectedHistory"],data["downstream_flows"]))
    # push!(plots, plotRestrictionsGrid(data["ρsHistory"]))

    # Load the image and add it to the plots
    
    # img_file = drawNetworkPNG(meta.populations,net.connections,infectedHistory, susceptibleHistory, ρsHistory)
    # img_file = drawNetworkKarnak(meta, data)
    # animate_network(meta.populations,convert(Matrix{Float64}, adjacency_matrix(meta.S.net.graph)),infectedHistory, susceptibleHistory, recoveredHistory, ρsHistory)
    # img = load(img_filename)
    # img_plot = plot(img, seriestype=:image)
    # push!(plots, (img_plot,800))
    if save
        savePlots(plots, meta)
    end
    return alignAndCombinePlots(plots)
end


# TODO create a function that unifies the functionality of all these functions
function plotInfEvolution(data; yLog=false, xLog = false)
    title = "Prevalence of Infected" * (xLog ? " (Log-x)" : "") * (yLog ? " (Log-y)" : "")
    xAxis = data["days"]
    infected_percent,avg_infected_percent = data["infectedHistory"],data["infectedAvgHistory"]
    nTimeSteps, nPopulations = size(infected_percent)

    infected_cleaned = replace(x -> x <= 0 ? NaN : x, infected_percent)
    colors = palette(:jet, nPopulations)
    infected_percent = yLog ? infected_cleaned : infected_percent
    color_index = 1
    p = plot()
    yLog && plot!(p, yscale=:log10)
    xLog && plot!(p, xscale=:log10)
    for i in 1:nPopulations
        plot!(p,xAxis, infected_percent[:, i], label="Pop. $i", color=colors[color_index])
        color_index += 1
    end

    plot!(p, xAxis, avg_infected_percent, label="Average", linestyle=:dash, linewidth=2)
    xlabel!(p,"Time (days)")
    ylabel!(p,"Infected Population Fraction")
    title!(p, title)
    height = yLog ? 600 : 400
    return p, height, title
end


function plotRestrictions(data)
    title = "Evolution of Path Mobility Restrictions"
    xAxis = data["days"]
    nTimeSteps, nρs= size(data["downstream_ρs"])
    p = plot(legend=false)
    
    colors = palette(:jet, nρs)
    
    color_index = 1
    for i in 1:nρs
            plot!(p, xAxis, data["downstream_ρs"][:,i] , label="P($i,$(i+1))",
             color=colors[color_index], ylim=(0,1))
            color_index += 1
    end
    
    xlabel!(p, "Time (days)")
    ylabel!(p, "Mobility Restrictions")
    title!(p, title)
    return p, 400, title
end

function plotInfectedFlow(data)
    title = "Evolution of Infected Flow"
    xAxis = data["days"]
    downstream_flows = data["downstream_flows"]
    nTimeSteps, nFlows = size(downstream_flows)
    p = plot(legend=false)
    
    colors = palette(:jet, nFlows)
    cleanDownstreamFlows = replace(x -> x <= 0 ? NaN : x, downstream_flows)

    color_index = 1
    for i in 1:nFlows
        plot!(p, xAxis, cleanDownstreamFlows[:,i], label="Flow($i,$(i+1))", color=colors[color_index],yscale=:log10)
        color_index += 1
    end
    
    xlabel!(p, "Time (days)")
    ylabel!(p, "Infected Flow")
    title!(p, title)
    return p, 400, title
end


function plotAvgRestrictions(ρsAvgHistory)
    title = "Evolution of Mobility Restrictions"
    
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
    title!(p, title)
    return p, 600, title
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
    
    return p, 800, ""
end

function plotInfectionDays(data,epi)
    title = "Spread and Peak Infection Days"

    pathLengths = data["pathLengths"]

    p = plot(xlabel="Path Length", ylabel="day", 
                legendfontsize=12, mc=:blue,width=800, size=(800,600))
    scatter!(p,pathLengths, data["spreadInfInd"], label="Spread Infection Day")
    epi.γ>0 && scatter!(p, pathLengths, data["peakInfInd"], label="Peak Infection Day", mc=:red)
    scatter!(p, pathLengths, data["firstOrderSol"], label="First Order Sol.", mc=:yellow)
    plot!(p, pathLengths, data["firstOrderApprox"], label="Lambert W rate", mc=:white)
    title!(p, title)
    return p, 800, title
end

function plotSpreadRates(data)
    pathLengths = data["pathLengths"]
    spreadInfInd = data["spreadInfInd"]

    # Calculate the inverse of the difference of the number of days between consecutive nodes
    inv_diff_days = 1.0 ./ (spreadInfInd[2:end] .- spreadInfInd[1:end-1])
    inv_diff_firstOrderSol = 1.0 ./ (data["firstOrderSol"][2:end] .- data["firstOrderSol"][1:end-1])
    inv_diff_firstOrderApprox = 1.0 ./ (data["firstOrderApprox"][2:end] .- data["firstOrderApprox"][1:end-1])
    
    # Create a scatter plot
    p = scatter(pathLengths, data["1/Δt"], label="1/Δt",
                xlabel="Path Length", ylabel="Inverse Δt", size=(800,600),
                legendfontsize=15, tickfontsize=15, guidefontsize=15, left_margin=10Plots.mm,
                mc=:blue)
    scatter!(p, pathLengths[2:end], inv_diff_firstOrderSol, label="First Order Sol.", mc=:yellow)
    plot!(p, pathLengths[2:end], inv_diff_firstOrderApprox, label="First Order Approx.")
    # make y-axis start from 0
    plot!(p, ylims=(0, maximum(inv_diff_days)*1.1))
    title = "Spread Rates"
    title!(p, title)
    return p, 600, title
end


function plot_cumulative_flow(data::Dict)
    cumuInfMob =[data["cumuInfMob"][i, i+1] for i in 1:size(data["cumuInfMob"], 1)-1]
    nPopulations = length(cumuInfMob)
    p = scatter(1:nPopulations-1, cumuInfMob, legend=false, yscale=:log10)
    plot!(p, xlims=(0, nPopulations))
    xlabel!(p, "Path length")
    ylabel!(p, "Cumulative Infected Flow")
    title= "Cumulative Infected Flow"
    title!(p, title)
    return p, 800, title
end

function plot_consecutive_infected(infectedHistory::Array{Float64,2}; log_scale::Bool=false)
    
    nTimeSteps, nPopulations = size(infectedHistory)
    p = plot()
    colors = palette(:jet, nPopulations)
    color_index = 1
    for i in 1:nPopulations-1
        phase_shift = (i - 1)/nPopulations  # Adjust phase shift as needed
        plot!(p, infectedHistory[2:end,i], infectedHistory[2:end,i+1], label = "($i,$(i+1))", color = colors[color_index], linestyle=:dash, phase=phase_shift)
        color_index += 1
    end
    xlabel!(p, "Population i infected prevalence")
    ylabel!(p, "Population i+1 infected prevalence")
    title = "Consecutive Infecteds prevalence" * (log_scale ? " (Log Scale)" : "")
    title!(p, title)
    
    if log_scale
        plot!(p,xscale=:ln,yscale=:ln)
    end
    
    return p, 800, title
end

function plot_infect_ρ(infectedHistory::Array{Float64,2}, ρsHistory::Array{Float64,3})
    
    nTimeSteps, nPopulations = size(infectedHistory)
    p = plot()
    colors = palette(:jet, nPopulations)
    color_index = 1
    cleanInfectedHistory = replace(x -> x <= 0 ? NaN : x,infectedHistory)
    cleanρsHistory = replace(x -> x <= 0 ? NaN : x, ρsHistory)
    for i in 1:nPopulations-1
        plot!(p, cleanInfectedHistory[2:end,i], cleanρsHistory[2:end,i+1,i], label = "($i,$(i+1))", color=colors[color_index], xscale=:log10)
        color_index += 1
    end
    plot!(p,legend=false)
    xlabel!(p, "Population i infected prevalence")
    ylabel!(p, "Population i+1 infected prevalence")
    title = "Evolution of infected prevalence and ρ"
    title!(p, title)
    return p, 800, title
end

function plot_infect_flow(infectedHistory, downstream_flows)
    nTimeSteps, nPopulations = size(infectedHistory)
    p = plot()
    colors = palette(:jet, nPopulations)
    cleanInfectedHistory = replace(x -> x <= 0 ? NaN : x, infectedHistory)
    cleanDownstreamFlows = replace(x -> x <= 0 ? NaN : x, downstream_flows)
    color_index = 1
    for i in 1:nPopulations-1
        plot!(p,  cleanDownstreamFlows[2:end, i],infectedHistory[2:end, i+1], label = "($i,$(i+1))", color=colors[color_index], xscale=:log10,yscale=:log10)
        color_index += 1
    end
    
    plot!(p, legend=false)
    ylabel!(p, "Population i infected prevalence")
    xlabel!(p, "Population i+1 infected flow rate")
    title = "Evolution of infected flow rate"
    title!(p, title)
    return p, 800, title
end