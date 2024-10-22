Base.@kwdef mutable struct SIRS_epidemic
    γ::Float64 # recovery rate
    β::Float64 # internal population infection rate
    σ::Float64 # immunity loss rate
    μ::Float64 # external link population infection rate
end


