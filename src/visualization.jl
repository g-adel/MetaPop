
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
    # angle2 = atan(p2.y - center.y, p2.x - center.x)

    return center, radius, angle1
end

function drawPopulation(pop::Population;radiusScale = 15)
    # Calculate the fractions
    maxRadius = sqrt(pop.size)*radiusScale

    Luxor.translate(pop.position)
    Luxor.rotate(-π/2)
    sAngle = 2π- pop.S* 2 * pi #cw
    iAngle = 2π - pop.I * 2 * pi

    setcolor("blue")
    sector(O, 0, maxRadius, 0, sAngle, :fill)
    setcolor("red")
    sector(O, 0, maxRadius, sAngle, sAngle + iAngle, :fill)
    setcolor("green")
    sector(O, 0, maxRadius, sAngle + iAngle, 2 * π, :fill)
        
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

function drawPopulationHistory(pop, popSusceptibleHistory, popInfectedHistory,)
    nTimeSteps = length(popSusceptibleHistory)
    popC = deepcopy(pop)
    nPts = 40
    Luxor.translate(pop.position+Point(-nPts/2,-nPts/2))
    for r in 1:nPts
        S = popSusceptibleHistory[floor(Int,r*nTimeSteps/nPts)]
        I = popInfectedHistory[floor(Int,r*nTimeSteps/nPts)]
        R = 1 - S - I
        setcolor(0,0,1)
        rect(Point(r,0),1.5,nPts*S,:fill)
        setcolor("green")
        rect(Point(r,nPts*(S)),1.5,nPts*R,:fill)
        setcolor("red")
        rect(Point(r,nPts*(S+R)),1.5,nPts*I,:fill)
        
    end
    origin()
end 

function drawConnection(i,j,populations; PNG = false)
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
    mobRestriction1 = populations[startInd].restrictions[endInd]
    mobRestriction2 = populations[endInd].restrictions[startInd]
    nSegs = 20
    center, radius, angle1 = drawArc2PntAngle(populations[startInd].position,populations[endInd].position, angle)

    
    # lineStep = (populations[j].position-populations[i].position)/nSegs
    angleStep = angle/nSegs
    for k in 1:nSegs
        x = (k-1)/nSegs
        redness = infFrac1*(1-x) + infFrac2*x
        # transparency = max(mobRate1*(1-x)^2, mobRate2*(x)^2) * abs(nSegs*.5-k)/nSegs*4
        transparency::Float64 = parabolaS0(x,mobRestriction1,mobRestriction2)
        if (PNG)
            setcolor(0,0,0)
        else
            setcolor(1-(1-redness)^1.1,0,0,transparency)
        end
        arc(center,radius,angle1+angleStep*(k-1),angle1+angleStep*k,:stroke)
        # line(populations[i].position+lineStep*(k-1), populations[i].position+lineStep*k,:stroke)
    end
end

function drawConnections(populations,connections;PNG = false)
    setline(2) # Set line width
    sethue("black") # Set line color
    nPopulations = length(populations)
    for i in 1:nPopulations
        for j in 1:i
            if connections[i,j] == 1
                drawConnection(i,j,populations; PNG = PNG)
            end
        end
    end
end

function drawNetwork(populations,connections)
    background("white")
    origin()
    drawConnections(populations,connections)
    for pop in populations
        drawPopulation(pop)
    end
    finish()
end

function frame(scene::Scene, framenumber::Int, populations,infectedHistory, susceptibleHistory, recoveredHistory, restrictionsHistory,connections)
    # locally available variables: populations,connections,infectedHistory, restrictionsHistory
    # create `populationsSnapshot` based on framenumber
    populationsSnapshot = deepcopy(populations)
    for population in populationsSnapshot
        population.I = infectedHistory[framenumber, population.index]
        population.S = susceptibleHistory[framenumber, population.index]
        population.R = recoveredHistory[framenumber, population.index]
        population.restrictions = restrictionsHistory[framenumber, population.index, :]
    end

    drawNetwork(populationsSnapshot,connections)
end


function drawNetworkPNG(populations,connections,infectedHistory, susceptibleHistory,restrictionsHistory;filename = "MetaPopNet.png")
    @png begin
        background("white")
        origin()
        drawConnections(populations,connections;PNG = true)
        for pop in populations
            drawPopulationHistory(pop, susceptibleHistory[:,pop.index], infectedHistory[:,pop.index])
        end
        finish()
    end 400 400 filename
    return filename
end

function animate_network(populations,connections,infectedHistory, susceptibleHistory, recoveredHistory, restrictionsHistory;filename = "preview.gif")
    nFrames = size(infectedHistory,1)
    mymovie=Movie(W,H,"test",1:nFrames)
    Luxor.animate(mymovie,[Scene(mymovie,(s, f) -> frame(s, f, populations,infectedHistory, susceptibleHistory, recoveredHistory, restrictionsHistory,connections),1:nFrames)],creategif=true, framerate=2,pathname=filename)
end