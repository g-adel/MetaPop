using Luxor, Revise
using Plots
include("structs.jl")
include("populations.jl")
include("networkViz.jl")
include("network.jl")
include("plots.jl")
include("strategies.jl")

α::Float64 = 0.05 # recovery rate
β::Float64 = 0.5 # internal population infection rate
μ::Float64 = 0.02 # external link population infection rate


k_bar::Int = 4
nPopulations::Int = 11
populations::Array{Population, 1} = Array{Population, 1}(undef, nPopulations)
connections::Array{Float64, 2} = zeros(Float64, nPopulations, nPopulations)

nTimeSteps::Int = 100

initializePopulations!(populations)
fillConnectionMatrix!(connections,nPopulations, k_bar)

susceptible = zeros(nTimeSteps, nPopulations)
infected = zeros(nTimeSteps, nPopulations)
recovered = zeros(nTimeSteps, nPopulations)

for i in 1:nTimeSteps
    updateNetwork!(populations, connections)
    for j in 1:nPopulations
        susceptible[i, j] = populations[j].S
        infected[i, j] = populations[j].I
        recovered[i, j] = populations[j].R
    end
end
draw_network(populations,connections)


plotTimeEvolution(infected)
