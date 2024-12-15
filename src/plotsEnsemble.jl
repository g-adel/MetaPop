function plotEnsemble(datas,Ss,meta;save=false)
    plots = []
    # # β and μ
    # push!(plots, plot_beta_mu_spread_rates(datas,Ss))
    # push!(plots, plot_beta_mu_spread_rates(datas,Ss;plot_t3=true))
    
    # λ
    push!(plots, plot_lambda_avg_spread_rate(datas,Ss))
    push!(plots, plot_lambda_path_spread_rates(datas,Ss))
    push!(plots, plot_lambda_F₂(datas,Ss))
    push!(plots, plot_lambda_t₂(datas,Ss))
    push!(plots, plot_lambda_t₃(datas,Ss))
    if save
        savePlots(plots, meta)
    end
    return alignAndCombinePlots(plots)
end


function plot_beta_mu_spread_rates(datas::Array{Dict, 2},Ss;plot_t3=false)
    nRows, nCols = size(datas)
    
    # Initialize arrays to store the values
    referenceSpreadRates = [[] for _ in 1:nCols]
    firstOrderSolSpreadRates = [[] for _ in 1:nCols]
    LambertApproxSpreadRates = [[] for _ in 1:nCols]
    Δt₃ = [[] for _ in 1:nCols]
    

    spreadRates = [[] for _ in 1:nCols]
    βs=[]
    μs=[]
    
    # Extract values from each dictionary in the 2D array
    for i in 1:nRows
        for j in 1:nCols
            push!(referenceSpreadRates[j], datas[i, j]["referenceSpreadRate"])
            push!(firstOrderSolSpreadRates[j], datas[i, j]["firstOrderSolSpreadRate"])
            push!(LambertApproxSpreadRates[j], 1 ./firstOrderApprox(Ss[i,j].epi))
            push!(spreadRates[j], datas[i, j]["avgSpreadRate"])
            push!(Δt₃[j], 1/(datas[i, j]["spreadInfInd"][3]-datas[i, j]["spreadInfInd"][2]))
            i==1 && push!(βs,Ss[1,j].epi.β)
        end
        push!(μs,Ss[i,1].epi.μ) 
    end
    @show βs μs

    min_μ = μs[1]
    max_μ = μs[end]
    μs_cont = range(min_μ, max_μ, length=100)
    epis = [[deepcopy(Ss[1,1].epi) for _ in eachindex(μs_cont)]]
    for j in 2:nCols
        push!(epis,[deepcopy(Ss[1,j].epi) for _ in eachindex(μs_cont)])
    end
    Lambert_cont = zeros(Float64,nCols,100)
    log_cont = zeros(Float64,nCols,100)
    
    for i in 1:100
        for j in 1:nCols
            epis[j][i].μ=μs_cont[i]
            Lambert_cont[j,i] = 1/firstOrderApprox(epis[j][i])
            κ = epis[j][i].β-2*epis[j][i].μ-epis[j][i].γ
            log_cont[j,i] = 1/(1/κ*log(1+κ/epis[j][i].μ))
        end
    end

    p = plot(xlabel="μ", ylabel="Spread Rate (nodes/day)", title="Spread Rates",legendfontsize=10)
    labels = ["Actual Rate","Equation Sol.","Lambert W","1/Δt₃","log"]
    for j in 1:nCols
        if j>1 labels=["","","","",""] end

        β = round(βs[j], digits=1)
        if !plot_t3
            plot!(p, μs, spreadRates[j], lw=2, lc=:green,linestyle=:dash,label="")
            scatter!(p, μs, spreadRates[j], label=labels[1], lw=2, mc=:green)
            # scatter!(p, μs, referenceSpreadRates[j], label="Ref Rate (β= $(βs[j]))", lw=2, mc=:black)
            plot!(p, μs, firstOrderSolSpreadRates[j],linestyle=:dash, lc=:red,label="")
            scatter!(p, μs, firstOrderSolSpreadRates[j], label=labels[2], lw=2, mc=:red)
            # scatter!(p, μs, firstOrderSolSpreadRates[j], label="My Rate (β= $(βs[j]))", lw=2, mc=:red)
            # plot!(p, μs, firstOrderApproxSpreadRates[j], label="Approx Rate (β= $(βs[j]))", lw=2, lc=:blue)
            annotate!(μs[end]*.8, spreadRates[j][end]*1, "β=$(β)")
        else
            plot!(p, μs, Δt₃[j,:], lw=1, lc=:green,linestyle=:dash,label="")
            scatter!(p, μs, Δt₃[j,:], mc=:green,label=labels[4])
            plot!(p, μs_cont,log_cont[j,:], label=labels[5], lw=2, lc=:red)
            annotate!(μs[end]*.8, Δt₃[j][end]*1, "β=$(β)")
        end
        plot!(p, μs_cont,Lambert_cont[j,:], label=labels[3], lw=2, lc=:blue)

    end
    
    title = plot_t3 ? "Δt₃ for different μ and β" : "Spread Rates for different μ and β"
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
    p = plot(xlabel="λ", ylabel="t₂",legendfontsize=15, legend=:topleft, xscale=:log10)
    scatter!(p, λs, t₂s_cleaned,label="t₂", mc=:yellow)
    plot!(p, λs,t₂Estimates_cleaned ,label="t₂ Estimate", lw=4)
    title = "t₂ vs λ (log-x)"
    title!(p, title)
    return p, 800, title
end

function plot_lambda_t₃(datas::Array{Dict, 1}, Ss)
    n = length(datas)
    
    λs = []
    t₃s = []
    t₃Estimates = []
    
    for i in 1:n
        S=Ss[i]
        push!(t₃s,datas[i]["spreadInfInd"][3])
        push!(λs, S.strat.λ)
        # push!(t₃Estimates, log(S.sim.I₀ * S.strat.λ)/(S.epi.β))
    end
    
    t₃s_cleaned = [t₃ <= 0 ? NaN : t₃ for t₃ in t₃s]
    t₃Estimates_cleaned = [t₃Estimates[i] <= 0 ? NaN : t₃Estimates[i] for i in 1:n]
    p = plot(xlabel="λ", ylabel="t₃",legendfontsize=15, legend=:topleft, xscale=:log10)
    scatter!(p, λs, t₃s_cleaned,label="t₃", mc=:yellow)
    # plot!(p, λs,t₃Estimates_cleaned ,label="t₃ Estimate", lw=4)
    title = "t₃ vs λ (log-x)"
    title!(p, title)
    return p, 800, title
end