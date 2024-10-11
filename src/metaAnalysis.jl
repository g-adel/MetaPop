function multiScenario(S)
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