# return array of floats of size nPopulations for mobility rate 

function linearUniformResponse(index,populations,localConnections::Array{Float64, 1})
   # compute average incoming and outgoing infection rate
    totalIncomingInfected = 0;    totalOutgoingInfected = 0; 
    nLinks = 0
    for (i, connectionStrength) in enumerate(localConnections)
        if connectionStrength > 0
            totalIncomingInfected += connectionStrength * μ * populations[i].I
            totalOutgoingInfected += connectionStrength * μ * populations[index].I
            nLinks+=1
        end
    end
    netRate = totalIncomingInfected - totalOutgoingInfected
    connectivity = 1-netRate/(nLinks)
    return ones(length(populations))*clamp(connectivity,0,1)
end
