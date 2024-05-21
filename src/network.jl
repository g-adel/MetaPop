
function fillConnectionMatrix!(connections::Array{Float64, 2}, nPopulations::Int, k_bar)
    for i in 1:nPopulations
        columns = [mod(i+j-1, nPopulations)+1 for j in -k_bar÷2:k_bar÷2 if mod(i+j-1, nPopulations)+1 != i]
        connections[i, columns] .= 1
    end
end


function updateNetwork!(populations, connections)
    populations_copy = deepcopy(populations)

    for (i, population_copy) in enumerate(populations_copy)
       updatePopulation!(population_copy, connections[i,:],populations_copy)

        populations[i] = population_copy
    end
end