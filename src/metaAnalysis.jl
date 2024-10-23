function multiScenario_μβ(S)
    nPtsX=nPtsY=5
    βs=range(0.1,1,nPtsX)
    μs=range(0.01,0.1,nPtsY)
    Ss=Array{Scenario,2}(undef,nPtsX,nPtsY)
    for i in 1:nPtsX
        for j in 1:nPtsY
            epi=deepcopy(S.epi)
            epi.μ=μs[i]
            epi.β=βs[j]
            Ss[i,j]=deepcopy(S)
            Ss[i,j].epi=epi
        end
    end
    return Ss
end

function multiScenario_λ(S)
    nPtsX=5
    λs=range(0.1,1,nPtsX)
    Ss=Array{Scenario,1}(undef,nPtsX)
    for i in 1:nPtsX
        strat=deepcopy(S.strat)
        strat.λ=λs[i]
        Ss[i]=deepcopy(S)
        Ss[i].strat=strat
    end
    return Ss
end