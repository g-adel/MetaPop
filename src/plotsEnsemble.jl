function plotMulti(datas,Ss,meta;save=false)
    plots = []
    push!(plots, plotInfectionDaysEnsemble(datas,meta.S.epi))
    push!(plots, plotInfectionDaysEnsemble(datas,meta.S.epi;distance="resistance_distances"))
    if save
        savePlots(plots, meta)
    end
    return alignAndCombinePlots(plots)
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

        scatter!(p,pathLengths, data["spreadInfInd"], label=labels[1], markerstrokecolor=:green,markercolor=:white)
        epi.γ>0 && scatter!(p, pathLengths, data["peakInfInd"], label=labels[2], markerstrokecolor	=:red,markercolor=:white)
        
        color_index+=1
        distance=="pathLengths" && plot!(p, pathLengths, data["firstOrderApprox"], label=labels[3], color=:blue)
        
        distance=="resistance_distances" && scatter!(p, pathLengths, data["firstOrderApproxResistance"],label=labels[3], mc=:blue)
        Labels=["","",""]
    end

    return p, 800, title
end

function plotInfectionDaysStrategies(datas,epi;distance="pathLengths")
    title = "Spread and Peak Infection Days ($distance)"
    p = plot( ylabel="Spread infection day", 
                legendfontsize=12, mc=:blue,width=800, size=(600,600))
    plot!(p,xlabel=distance)
    colors = palette(:jet, length(datas))
    color_index = 1
    labels=["Spread Infection Day","Peak Infection Day", "Lambert W rate"]
    for data in datas
        
        pathLengths = data[distance]

        scatter!(p,pathLengths, data["spreadInfInd"], label=labels[1], markerstrokecolor=:green,markercolor=:white)
        epi.γ>0 && scatter!(p, pathLengths, data["peakInfInd"], label=labels[2], markerstrokecolor	=:red,markercolor=:white)
        
        color_index+=1
        distance=="pathLengths" && plot!(p, pathLengths, data["firstOrderApprox"], label=labels[3], color=:blue)
        
        distance=="resistance_distances" && scatter!(p, pathLengths, data["firstOrderApproxResistance"],label=labels[3], mc=:blue)
        Labels=["","",""]
    end

    return p, 800, title
end