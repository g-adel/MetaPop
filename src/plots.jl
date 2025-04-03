function alignAndCombinePlots(plots)
    gr()
    height = sum(p[2] for p in plots)
    for p in plots
        plot!(p[1], left_margin=20Plots.mm, bottom_margin=10Plots.mm, tickfontsize=10, guidefontsize=15)
        title!(p[1],p[3])
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
    push!(plots, plotInfEvolutionCumu(data))
    push!(plots, plotInfEvolution(data;yLog=true))
    if adaptive
        push!(plots, plotRestrictions(data))
        # push!(plots, plot_infect_ρ(data["infectedHistory"],data["ρsHistory"]))
        # push!(plots, plotAvgRestrictions(data["ρsAvgHistory"]))    
    end
    S.sim.critRange>0 && push!(plots, plotInfEvolution(data;yLog=true,xLog=true))
    # push!(plots, plotInfectedFlow(data))
    push!(plots, plotInfectionDays(data,S.epi))
    push!(plots, plotSpreadRates(data))
    push!(plots, plotFlows(data))
    push!(plots, plotFlowsCumu(data))
    # push!(plots, plot_cumulative_flow(data))
    # push!(plots, plot_consecutive_infected(data["infectedHistory"],log_scale=false))
    # (S.net.topology != PathGraph) && push!(plots, plot_image(drawNetworkKarnak(meta, data)))
    # push!(plots, plot_infect_flow(data["infectedHistory"],data["downstream_flows"]))
    # push!(plots, plotRestrictionsGrid(data["ρsHistory"]))
    
    # animate_network(meta.populations,convert(Matrix{Float64}, adjacency_matrix(meta.S.net.graph))
    # ,data["infectedHistory"], data["susceptibleHistory"], data["recoveredHistory"], data["ρsHistory"])
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
    yLog && plot!(p, yscale=:log10);    xLog && plot!(p, xscale=:log10)
    nPopulations>10 && plot!(p,legend=false)
    for i in 1:nPopulations
        plot!(p,xAxis, infected_percent[:, i], label="Pop. $i", color=colors[color_index])
        color_index += 1
    end

    plot!(p, xAxis, avg_infected_percent, label="Average", linestyle=:dash, linewidth=2, lc=:black)
    xlabel!(p,"Time (days)")
    ylabel!(p,"Infected Population Fraction")
    height = yLog ? 600 : 400
    return p, height, title
end

function plotInfEvolutionCumu(data; yLog=false, xLog = false)
    title = "Prevalence of Infected" * (xLog ? " (Log-x)" : "") * (yLog ? " (Log-y)" : "")
    xAxis = data["days"]
    infected_percent = data["infectedHistory"]
    nTimeSteps, nPopulations = size(infected_percent)
    infected_percent_cumu = zeros(nTimeSteps)

    p = plot(size=(600,600))
    infected_cleaned = replace(x -> x <= 0 ? NaN : x, infected_percent)
    colors = palette(:jet, nPopulations)
    infected_percent = yLog ? infected_cleaned : infected_percent
    color_index = 1
    yLog && plot!(p, yscale=:log10);    xLog && plot!(p, xscale=:log10)
    nPopulations>10 && plot!(p,legend=false)
    
    for i in 1:nPopulations
        plot!(p, xAxis, infected_percent_cumu, fillrange = (infected_percent_cumu.+infected_percent[:, i]), fillalpha = 0.35,label="")
        infected_percent_cumu.+=infected_percent[:, i]
        plot!(p,xAxis, infected_percent_cumu, label="", linewidth=0.5)
        color_index += 1
    end

    plot!(p, xAxis, infected_percent_cumu, label="Total", linestyle=:dash, linewidth=2, lc=:black)
    xlabel!(p,"Time (days)")
    ylabel!(p,"Infected Population Fraction")
    height = yLog ? 600 : 400
    return p, height, title
end


function plotRestrictions(data)
    title = "Evolution of Path Mobility Restrictions"
    xAxis = data["days"]
    nTimeSteps, nρs= size(data["downstream_ρs"])
    # p = plot(legend=false)
    p=plot(legend=:bottomright)
    colors = palette(:jet, nρs)
    
    color_index = 1
    for i in 1:nρs
            plot!(p, xAxis, data["downstream_ρs"][:,i] , label="(ρ($i,$(i+1)))",
             color=colors[color_index], ylim=(0,1))
            color_index += 1
    end
    
    xlabel!(p, "Time (days)")
    ylabel!(p, "Mobility Restrictions")
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
    return p, 400, title
end

function plotFlows(data)
    title = "Evolution of Mobility Flow"
    xAxis = data["days"][2:end]
    flowsHistory = data["flowsHistory"]
    nTimeSteps, nFlows = size(flowsHistory)
    totalFlow = zeros(nTimeSteps-1)
    p = plot(legend=false)
    # p=plot()
    colors = palette(:jet, nFlows)

    color_index = 1
    for i in 1:nFlows
        plot!(p, xAxis, flowsHistory[2:end,i], label="M_$i)", linewidth=0.5)
        totalFlow.+=flowsHistory[2:end,i]
        color_index += 1
    end
    plot!(p, xAxis, totalFlow./nFlows, label="Average", linestyle=:dash, linewidth=3, lc=:black)

    xlabel!(p, "Time (days)")
    ylabel!(p, "Mobility Flow")
    return p, 400, title
end

function plotFlowsCumu(data)
    title = "Cumulative Mobility Flow"
    xAxis = data["days"][2:end]
    flowsHistory = data["flowsHistory"]
    nTimeSteps, nFlows = size(flowsHistory)
    cumulativeFlow = zeros(nTimeSteps-1)

    p = plot(size=(600, 600), legend=false)
    colors = palette(:jet, nFlows)
    color_index = 1

    for i in 1:nFlows
        plot!(p, xAxis, cumulativeFlow, fillrange=(cumulativeFlow .+ flowsHistory[2:end, i]), fillalpha=0.35, label="")
        cumulativeFlow .+= flowsHistory[2:end, i]
        plot!(p, xAxis, cumulativeFlow, label="", linewidth=0.5)
        color_index += 1
    end

    plot!(p, xAxis, cumulativeFlow, label="Total", linestyle=:dash, linewidth=2, lc=:black)
    xlabel!(p, "Time (days)")
    ylabel!(p, "Cumulative Mobility Flow")
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
    return p, 600, title
end


function plot_cumulative_flow(data::Dict)
    title= "Cumulative Infected Flow"
    cumuInfMob =[data["cumuInfMob"][i, i+1] for i in 1:size(data["cumuInfMob"], 1)-1]
    nPopulations = length(cumuInfMob)
    p = scatter(1:nPopulations-1, cumuInfMob, legend=false, yscale=:log10)
    plot!(p, xlims=(0, nPopulations))
    xlabel!(p, "Path length")
    ylabel!(p, "Cumulative Infected Flow")
    
    return p, 800, title
end

function plot_consecutive_infected(infectedHistory::Array{Float64,2}; log_scale::Bool=false)
    title = "Consecutive Infecteds prevalence" * (log_scale ? " (Log Scale)" : "")
    
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

    
    if log_scale
        plot!(p,xscale=:ln,yscale=:ln)
    end
    
    return p, 800, title
end

function plot_infect_ρ(infectedHistory::Array{Float64,2}, ρsHistory::Array{Float64,3})
    title = "Evolution of infected prevalence and ρ"
    
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
    return p, 800, title
end

function plot_infect_flow(infectedHistory, downstream_flows)
    title = "Evolution of infected flow rate"
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
    return p, 800, title
end

function plot_image(img_filename)
    p=plot(load(img_filename), seriestype=:image, framestyle=:none,size=(1000,1000))
    return p, 800, "Network"
end