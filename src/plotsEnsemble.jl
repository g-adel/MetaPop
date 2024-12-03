function plotEnsemble(datas,Ss,meta;save=false)
    plots = []
    push!(plots, plot_lambda_avg_spread_rate(datas,Ss))
    push!(plots, plot_lambda_path_spread_rates(datas,Ss))
    push!(plots, plot_lambda_F₂(datas,Ss))
    push!(plots, plot_lambda_t₂(datas,Ss))
    # combinedPlots = 
    if save
        savePlots(plots, meta)
    end
    return alignAndCombinePlots(plots)
end


function plot_beta_mu_spread_rates(datas::Array{Dict, 2},Ss)
    nRows, nCols = size(datas)
    
    # Initialize arrays to store the values
    referenceSpreadRates = [[] for _ in 1:nCols]
    firstOrderSolSpreadRates = [[] for _ in 1:nCols]
    firstOrderApproxSpreadRates = [[] for _ in 1:nCols]
    spreadRates = [[] for _ in 1:nCols]
    βs=[]
    μs=[]
    
    # Extract values from each dictionary in the 2D array
    for i in 1:nRows
        for j in 1:nCols
            push!(referenceSpreadRates[j], datas[i, j]["referenceSpreadRate"])
            push!(firstOrderSolSpreadRates[j], datas[i, j]["firstOrderSolSpreadRate"])
            push!(firstOrderApproxSpreadRates[j], 1 ./firstOrderApprox(Ss[i,j].epi))
            push!(spreadRates[j], datas[i, j]["avgSpreadRate"])
            i==1 && push!(βs,Ss[1,j].epi.β)
        end
        push!(μs,Ss[i,1].epi.μ) 
    end
    @show βs μs
    # Create the plot
    p = plot(xlabel="μ", ylabel="Rate (nodes/day)", title="Spread Rates")
    
    for j in 1:nCols
        plot!(p, μs, spreadRates[j], label="Actual Rate (β= $(βs[j]))", lw=2, lc=:green, legendfontsize=10)
        # scatter!(p, μs, referenceSpreadRates[j], label="Ref Rate (β= $(βs[j]))", lw=2, mc=:black)
        plot!(p, μs, firstOrderSolSpreadRates[j], label="My Rate (β= $(βs[j]))", lw=2, lc=:red)
        # scatter!(p, μs, firstOrderSolSpreadRates[j], label="My Rate (β= $(βs[j]))", lw=2, mc=:red)
        # plot!(p, μs, firstOrderApproxSpreadRates[j], label="Approx Rate (β= $(βs[j]))", lw=2, lc=:blue)
        scatter!(p, μs, firstOrderApproxSpreadRates[j], label="Approx Rate (β= $(βs[j]))", lw=2, mc=:blue)
    end
    title = "Spread Rates"
    title!(p, title)

    return p, 800, title
end

function plot_lambda_avg_spread_rate(datas::Array{Dict, 1}, Ss)
    n = length(datas)
    
    asympSpreadRates = []
    initSpreadRates = []
    avgSpreadRates = []
    λs = []
    
    for i in 1:n
        push!(asympSpreadRates, datas[i]["asympSpreadRate"])
        push!(initSpreadRates, datas[i]["initSpreadRate"])
        push!(avgSpreadRates, datas[i]["avgSpreadRate"])
        push!(λs, Ss[i].strat.λ)
    end
    
    # @show λs 
    p = plot(xlabel="λ", ylabel="Rate (nodes/day)", title="Spread Rates", ylims=(0,.2), xscale=:log10, legendfontsize=10)
    plot!(p, λs, asympSpreadRates, lw=1, lc=:green, linestyle=:dash, label="")
    scatter!(p, λs, asympSpreadRates, label="Asymptotic Spread Rate",lw=2, mc=:green)
    plot!(p, λs, initSpreadRates, lw=1, lc=:blue, linestyle=:dash, label="")
    scatter!(p, λs, initSpreadRates, label="Initial Spread Rate", lw=2, mc=:blue)
    plot!(p, λs, avgSpreadRates, lw=1, lc=:red, linestyle=:dash, label="")
    scatter!(p, λs, avgSpreadRates, label="Average Spread Rate", lw=2, mc=:red)

    # plot!(p, λs, firstOrderSolSpreadRates, label="My Rate", lw=2, lc=:red)
    # scatter!(p, λs, firstOrderApproxSpreadRates, label="Approx Rate", lw=2, mc=:blue)
    title = "Spread Rates"
    title!(p, title)

    return p, 800, title
end

function plot_lambda_path_spread_rates(datas::Array{Dict, 1}, Ss)
    n = length(datas)
    
    λs = []
    spreadRates = []
    
    # Extract values from each dictionary in the array
    for i in 1:n
        push!(spreadRates,datas[i]["1/Δt"])
        push!(λs, Ss[i].strat.λ)
    end
    
    pathLengths = datas[1]["pathLengths"]

    colors = palette(:jet, n)
    color_index = 1

    p= scatter(legendfontsize=10,xlabel="Path Length", ylabel="Inverse Δt", size=(1000,800))
    for i in 1:n
        plot!(p,pathLengths[2:end], spreadRates[i][2:end], lw=1,lc=colors[color_index], linestyle=:dash, label="")
        scatter!(p, pathLengths[2:end], spreadRates[i][2:end], label="λ=$(λs[i])", mc=colors[color_index])
        color_index += 1
    end
    scatter!(p, ylims=(0, .2))
    title =  "Inverse Δt vs Path Length"
    title!(p, title)
    return p, 800, title
end


function plot_lambda_F₂(datas::Array{Dict, 1}, Ss)
    n = length(datas)
    
    λs = []
    F₂s = []
    cs=[]
    
    for i in 1:n
        S=Ss[i]
        push!(F₂s,datas[i]["infectedHistory"][2,2])
        push!(λs, S.strat.λ)
        # push!(F₂Estimates,1/λs[i])
        cs=push!(cs,S.strat.λ*S.epi.μ*S.sim.I₀/S.epi.β)
    end
    F₂Estimates = 1 ./ λs
    # F₂_low_c = 1 ./ (λs .+ 1 ./ cs)
    # @show F₂_low_c
    p = plot(legendfontsize=15, xscale=:log10, yscale=:log10)
    scatter!(p, λs, F₂s,label="F₂", mc=:yellow)
    plot!(p, xlabel="λ", ylabel="F₂")
    plot!(p, λs,F₂Estimates,label="F₂ Estimate", lw=4)

    title = "F₂ vs λ (log-log)"
    title!(p, title)

    return p, 800, title
end

function plot_lambda_t₂(datas::Array{Dict, 1}, Ss)
    n = length(datas)
    
    λs = []
    t₂s = []
    t₂Estimates = []
    
    for i in 1:n
        S=Ss[i]
        push!(t₂s,datas[i]["spreadInfInd"][2])
        push!(λs, S.strat.λ)
        push!(t₂Estimates, log(S.sim.I₀ * S.strat.λ)/(S.epi.β))
    end
    @show t₂Estimates t₂s
    t₂s_cleaned = [t₂ <= 0 ? NaN : t₂ for t₂ in t₂s]
    t₂Estimates_cleaned = [t₂Estimates[i] <= 0 ? NaN : t₂Estimates[i] for i in 1:n]
    p = plot()
    scatter!(p, λs, t₂s_cleaned,label="t₂", mc=:yellow)
    plot!(p,xlabel="λ", ylabel="t₂",legendfontsize=15, legend=:topleft, xscale=:log10)
    plot!(p, λs,t₂Estimates_cleaned ,label="t₂ Estimate", lw=4)
    title = "t₂ vs λ (log-x)"
    title!(p, title)
    return p, 800, title
end