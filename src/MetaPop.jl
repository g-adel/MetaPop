using Luxor, Images, Revise
using Plots
foreach(include, ["util.jl","structs.jl", "epidemics.jl", 
        "populations.jl", "networkViz.jl", "network.jl", "plots.jl", 
        "strategies.jl"])

epi = SIR_epidemic(0.03, 0.2,0.01)
net = Network(11,4)
populations::Array{Population, 1} = Array{Population, 1}(undef, net.nPopulations)
connections::Array{Float64, 2} = zeros(Float64, net.nPopulations, net.nPopulations)

nTimeSteps::Int = 200

initializePopulations!(populations)
fillConnectionMatrix!(connections,net.nPopulations, net.k_bar)
susceptibleHistory = zeros(nTimeSteps, net.nPopulations)
infectedHistory = zeros(nTimeSteps, net.nPopulations)
recoveredHistory = zeros(nTimeSteps, net.nPopulations)
mobilityRatesHistory = zeros(nTimeSteps, net.nPopulations, net.nPopulations)

for i in 1:nTimeSteps
    updateNetwork!(populations, connections,epi)
    for j in 1:net.nPopulations
        susceptibleHistory[i, j] =  populations[j].S
        infectedHistory[i, j] = populations[j].I
        recoveredHistory[i, j] = populations[j].R
        mobilityRatesHistory[i, j, :] = populations[j].mobilityRates
    end
end


# draw_network(populations,connections)
# animate_network(populations,connections,infectedHistory, mobilityRatesHistory)

p1 = plotTimeEvolution(infectedHistory)
p2 = plotTotalConnectivity(mobilityRatesHistory,connections)
plot(p1, p2, layout = (2, 1), size = (500, 800))