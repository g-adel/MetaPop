
W = 400; H= 400;

function drawArc2PntAngle(p1::Point, p2::Point, angle::Float64)
    midpoint = 0.5 * (p1 + p2)

    dist = mag(p2 - p1)

    radius = dist / (2 * sin(angle / 2))

    dist_to_center = sqrt(radius^2 - (dist / 2)^2)

    direction = (p2 - p1) / dist
    direction = Point(-direction.y, direction.x)

    center = midpoint + dist_to_center * direction

    angle1 = atan(p1.y - center.y, p1.x - center.x)
    angle2 = atan(p2.y - center.y, p2.x - center.x)

    # arc(center, radius, angle1, angle2,:stroke)
    return center, radius, angle1
end

function drawPopulation(pop::Population)
    # Calculate the fractions
    total = pop.S + pop.I + pop.R
    s_fraction = pop.S / total
    i_fraction = pop.I / total
    r_fraction = pop.R / total
    maxRadius = sqrt(pop.size)*15

    Luxor.translate(pop.position)
    Luxor.rotate(-π/2)
    s_angle = s_fraction * 2 * pi
    i_angle = i_fraction * 2 * pi

    setcolor("blue")
    sector(O, 0, maxRadius, 0, s_angle, :fill)
    setcolor("red")
    sector(O, 0, maxRadius, s_angle, s_angle + i_angle, :fill)
    setcolor("green")
    sector(O, 0, maxRadius, s_angle + i_angle, 2 * pi, :fill)
        
    # function drawDisk(radius, color)
    #     setcolor(color)
    #     frac_radius = maxRadius * sqrt(radius)
    #     circle(O, frac_radius, :fill)
    # end
    
    # fractions = [(s_fraction, "green"), (i_fraction, "red"), (r_fraction, "blue")]
    # cumuFrac = 1
    # for (fraction, color) in fractions
    #     drawDisk(cumuFrac, color)
    #     cumuFrac -= fraction
    # end

    origin()
end

function drawConnection(i,j,populations)

    
    if (j>i)
        startInd, endInd = (j,i)
    else
        startInd, endInd = (i,j)
    end

    nPopulations = length(populations)
    angle = (abs(j-i))*π/(nPopulations/2)
    if (angle>π)
        angle = 2π-angle
        # print(angle)
        startInd, endInd = (endInd,startInd)
    end
    infFrac1 = populations[startInd].I
    infFrac2 = populations[endInd].I  
    mobRate1 = populations[startInd].mobilityRates[j]
    mobRate2 = populations[endInd].mobilityRates[i]
    nSegs = 20
    center, radius, angle1 = drawArc2PntAngle(populations[startInd].position,populations[endInd].position, angle)
    
    lineStep = (populations[j].position-populations[i].position)/nSegs
    angleStep = angle/nSegs
    for k in 1:nSegs
        redness::Float64 = infFrac1*(nSegs-k)/nSegs + infFrac2*(k-1)/nSegs
        transparency::Float64 = max(mobRate1*(nSegs-k)/nSegs, mobRate2*(k-1)/nSegs) * 2
        setcolor(redness^.5,0,0,transparency)
        arc(center,radius,angle1+angleStep*(k-1),angle1+angleStep*k,:stroke)
        # line(populations[i].position+lineStep*(k-1), populations[i].position+lineStep*k,:stroke)
    end
end

function drawConnections(populations,connections)
    setline(2) # Set line width
    sethue("black") # Set line color
    nPopulations = length(populations)
    for i in 1:nPopulations
        for j in 1:i
            if connections[i,j] == 1
                drawConnection(i,j,populations)
            end
        end
    end
end

function draw_network(populations,connections)
    # @png begin
        background("white")
        origin()
        drawConnections(populations,connections)
        for pop in populations
            drawPopulation(pop)
        end
        finish()
    # end 400 400 "MetaPopNet.png"
end

function frame(scene::Scene, framenumber::Int)
    # locally available variables: populations,connections,infectedHistory, mobilityRatesHistory
    # create `populationsSnapshot` based on framenumber
    populationsSnapshot = deepcopy(populations)
    for population in populationsSnapshot
        population.I = infectedHistory[framenumber, population.index]
        population.S = susceptibleHistory[framenumber, population.index]
        population.R = recoveredHistory[framenumber, population.index]
        population.mobilityRates = mobilityRatesHistory[framenumber, population.index, :]
    end

    draw_network(populationsSnapshot,connections)
end



function animate_network(populations,connections,infectedHistory, mobilityRatesHistory;filename = "preview.gif")
    nFrames = size(infectedHistory,1)
    mymovie=Movie(W,H,"test",1:nFrames)
    Luxor.animate(mymovie,[Scene(mymovie,frame,1:nFrames)],creategif=true, framerate=5,pathname=filename)
end