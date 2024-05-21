
W = 400; H= 400;

function drawPopulation(pop::Population)
    # Calculate the fractions
    total = pop.S + pop.I + pop.R
    s_fraction = pop.S / total
    i_fraction = pop.I / total
    r_fraction = pop.R / total
    maxRadius = sqrt(pop.size)*15

    
    function drawDisk(radius, color)
        setcolor(color)
        frac_radius = maxRadius * sqrt(radius)
        circle(O, frac_radius, :fill)
    end
    
    Luxor.translate(pop.position)
    
    fractions = [(s_fraction, "green"), (i_fraction, "red"), (r_fraction, "blue")]
    cumuFrac = 1
    for (fraction, color) in fractions
        drawDisk(cumuFrac, color)
        cumuFrac -= fraction
    end

    origin()
end

function drawConnection(i,j,populations)
    infFrac1 = populations[i].I
    infFrac2 = populations[j].I  
    minInf = min(infFrac1,infFrac2)
    # infFrac1 -= minInf; infFrac2 -= minInf;
    nSegs = 5
    lineStep = (populations[j].position-populations[i].position)/nSegs
    for k in 1:nSegs
        redness::Float64 = infFrac1*(nSegs-k)/nSegs + infFrac2*(k-1)/nSegs
        print(redness," - ")
        setcolor(redness*10,0,0)

        line(populations[i].position+lineStep*(k-1), populations[i].position+lineStep*k,:stroke)
    end
    strokepath()
end

function drawConnections(populations,connections)
    setline(2) # Set line width
    sethue("black") # Set line color
    for i in 1:nPopulations
        for j in i:nPopulations
            if connections[i,j] == 1
                drawConnection(i,j,populations)
            end
        end
    end
end

function draw_network(populations,connections)
    @png begin
        origin()
        drawConnections(populations,connections)
        for pop in populations
            drawPopulation(pop)
        end
        finish()
    end 400 400 "MetaPopNet.png"
end
