Base.@kwdef struct SIR_epidemic
    γ::Float64 # recovery rate
    β::Float64 # internal population infection rate
    μ::Float64 # external link population infection rate
end

Base.@kwdef struct SIRS_epidemic
    γ::Float64 # recovery rate
    β::Float64 # internal population infection rate
    σ::Float64 # immunity loss rate
    μ::Float64 # external link population infection rate
end

struct SEIR_epidemic
    γ::Float64 # recovery rate
    β::Float64 # internal population infection rate
    ε::Float64 # exposed to infected rate
    μ::Float64 # external link population infection rate
end

struct SEIRS_epidemic
    γ::Float64 # recovery rate
    β::Float64 # internal population infection rate
    ε::Float64 # exposed to infected rate
    σ::Float64 # immunity loss rate
    μ::Float64 # external link population infection rate
end





