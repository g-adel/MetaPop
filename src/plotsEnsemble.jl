function plotEnsemble(datas,Ss)
    plots = []
    push!(plots, plot_lambda_avg_spread_rate(datas,Ss))
    push!(plots, plot_lambda_path_spread_rates(datas,Ss))
    push!(plots, plot_lambda_F_2(datas,Ss))
    push!(plots, plot_lambda_t_2(datas,Ss))

    return combinePlots(plots)
end


function plot_multi_spread_rates(datas::Array{Dict, 2},Ss)
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
        plot!(p, μs, spreadRates[j], label="Actual Rate (β= $(βs[j]))", lw=2, lc=:green, legendfontsize=6)
        # scatter!(p, μs, referenceSpreadRates[j], label="Ref Rate (β= $(βs[j]))", lw=2, mc=:black)
        plot!(p, μs, firstOrderSolSpreadRates[j], label="My Rate (β= $(βs[j]))", lw=2, lc=:red)
        # scatter!(p, μs, firstOrderSolSpreadRates[j], label="My Rate (β= $(βs[j]))", lw=2, mc=:red)
        # plot!(p, μs, firstOrderApproxSpreadRates[j], label="Approx Rate (β= $(βs[j]))", lw=2, lc=:blue)
        scatter!(p, μs, firstOrderApproxSpreadRates[j], label="Approx Rate (β= $(βs[j]))", lw=2, mc=:blue)
    end
    
    return p, 800
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
    
    @show λs 
    # Create the plot
    p = plot(xlabel="λ", ylabel="Rate (nodes/day)", title="Spread Rates", ylims=(0,.2), xscale=:log10)
    plot!(p, λs, asympSpreadRates, label="Asymptotic Spread Rate", lw=2, lc=:green, legendfontsize=6)
    plot!(p, λs, initSpreadRates, label="Initial Spread Rate", lw=2, lc=:blue, legendfontsize=6)
    plot!(p, λs, avgSpreadRates, label="Average Spread Rate", lw=2, lc=:red, legendfontsize=6)

    # plot!(p, λs, firstOrderSolSpreadRates, label="My Rate", lw=2, lc=:red)
    # scatter!(p, λs, firstOrderApproxSpreadRates, label="Approx Rate", lw=2, mc=:blue)
    
    return p, 800
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

    p=plot(xlabel="Path Length", ylabel="Inverse Δt", size=(1000,800),
    legendfontsize=15, tickfontsize=15, guidefontsize=15)
    for i in 1:n
        plot!(p, pathLengths[2:end], spreadRates[i][2:end], label="λ=$(λs[i])", mc=:yellow)
    end
    plot!(p, ylims=(0, .2))

    return p, 800
end


function plot_lambda_F_2(datas::Array{Dict, 1}, Ss)
    n = length(datas)
    
    λs = []
    F_2s = []
    cs=[]
    
    for i in 1:n
        S=Ss[i]
        push!(F_2s,datas[i]["infectedHistory"][2,2])
        push!(λs, S.strat.λ)
        # push!(F_2Estimates,1/λs[i])
        cs=push!(cs,S.strat.λ*S.epi.μ*S.sim.I₀/S.epi.β)
    end
    F_2Estimates = 1 ./ λs
    F_2_low_c = 1 ./ (λs .+ 1 ./ cs)
    @show F_2_low_c
    p=plot(xlabel="λ", ylabel="F_2",legendfontsize=15, xscale=:log10, yscale=:log10)
    plot!(p, λs,F_2Estimates,label="F_2 Estimate", lw=4)
    scatter!(p, λs, F_2s,label="F_2", mc=:yellow)
    title!(p, "F_2 vs λ (log-log) ")
    return p, 800
end

function plot_lambda_t_2(datas::Array{Dict, 1}, Ss)
    n = length(datas)
    
    λs = []
    t_2s = []
    t_2Estimates = []
    
    for i in 1:n
        S=Ss[i]
        push!(t_2s,datas[i]["spreadInfInd"][2])
        push!(λs, S.strat.λ)
        push!(t_2Estimates, log(S.sim.I₀ * S.strat.λ)/(S.epi.β))
    end
    @show t_2Estimates t_2s
    t_2s_cleaned = [t_2 <= 0 ? NaN : t_2 for t_2 in t_2s]
    t_2Estimates_cleaned = [t_2Estimates[i] <= 0 ? NaN : t_2Estimates[i] for i in 1:n]

    p=plot(xlabel="λ", ylabel="t_2",legendfontsize=15, legend=:topleft, xscale=:log10)
    plot!(p, λs,t_2Estimates_cleaned ,label="t_2 Estimate", lw=4)
    scatter!(p, λs, t_2s_cleaned,label="t_2", mc=:yellow)
    title!(p, "F_2 vs λ (log-x) ")
    return p, 800
end