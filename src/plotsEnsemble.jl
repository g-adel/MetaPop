function plotEnsemble(datas,Ss,meta;save=false)
    plots = []
    # push!(plots, plotInfectionDaysStrategies(datas))
    # push!(plots, plotInfectionDaysStrategies(datas))
    # push!(plots, plotInfectionDaysEnsemble(datas,meta.S.epi))
    # push!(plots, plotInfectionDaysEnsemble(datas,meta.S.epi;distance="resistance_distances"))
    # push!(plots, plotInfectionDaysStrategies(datas))
    push!(plots, plotInfEvolutionEnsemble(datas))
    push!(plots, plotInfEvolutionEnsemble(datas))
    if save
        savePlots(plots, meta)
    end
    return alignAndCombinePlots(plots)
end

function plotInfEvolutionEnsemble(datas)
    title = "Prevalence of Infected (Ensemble)"
    
    p = plot(legend=false,size=(600,600))
    for data in datas
        xAxis = data["days"]
        avg_infected_percent = data["infectedAvgHistory"]
        plot!(p,xAxis, avg_infected_percent, label="")
    end
    xlabel!(p,"Time (days)")
    ylabel!(p,"Infected Population Size")
    return p, 600, title
end



function plotInfectionDaysEnsemble(datas,epi;distance="pathLengths")
    title = "Spread and Peak Infection Days ($distance)"
    p = plot( ylabel="Spread infection day", 
                legendfontsize=12, mc=:blue,width=800, size=(600,600))
    plot!(p,xlabel=distance)
    colors = palette(:jet, length(datas))
    color_index = 1
    labels=["Spread Infection Day","Peak Infection Day", "Lambert W rate"]
    for data in datas
        pathLengths = data[distance]

        scatter!(p,pathLengths, data["spreadInfInd"], label=labels[1], markercolor=:green,markerstrokecolor=colors[color_index])
        epi.γ>0 && scatter!(p, pathLengths, data["peakInfInd"], label=labels[2], markercolor	=:red,markerstrokecolor=colors[color_index])
        
        color_index+=1
        distance=="pathLengths" && plot!(p, pathLengths, data["firstOrderApprox"], label=labels[3], color=:blue)
        
        distance=="resistance_distances" && scatter!(p, pathLengths, data["firstOrderApproxResistance"],label=labels[3], mc=:blue)
        labels=["","",""]
    end

    return p, 800, title
end

function plotInfectionDaysStrategies(datas)
    title = "Spread Infection Days different strategies"
    p = plot(ylabel="Spread infection day", 
                legendfontsize=11, size=(700,500))
    plot!(p,xlabel="Path Length")
    colors = palette(:jet, length(datas))
    color_index = 1
    labels=["No Restriction","Global Proportional Restriction","Uniform Proportional Restriction",
    "Connection Proportional Restriction","Connection Logarithm Restriction"]
    i=1
    for data in datas
        strat = data["S"].strat.strategy
        pathLengths = data["pathLengths"]

        scatter!(p,pathLengths, data["spreadInfInd"], label=labels[i], markercolor=colors[color_index])
        plot!(p,pathLengths, data["spreadInfInd"], label="", color=colors[color_index],style=:dash)
        
        color_index+=1
 
        i+=1
    end

    return p, 700, title
end