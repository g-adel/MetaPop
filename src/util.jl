function structVals(s)
    return Tuple(getfield(s, name) for name in fieldnames(typeof(s)))
end

function parabolaS0(x,y1,y2)
    if(y2>y1)
        return y1+(x^2)*(y2-y1)
    else
        return y2+((1-x)^2)*(y1-y2)
    end
end

function mag(p::Point)
    return sqrt(p.x^2+p.y^2)
end

function perpendicular(p::Point)
    return Point(-p.y,p.x)
end

function perp(p::Point)
    return Point(-p.y,p.x)
end

function map(x,x0,x1,y0,y1)
    return (x-x0)/(x1-x0)*(y1-y0)+y0
end

function mapto01(x,x0,x1)
    return (x-x0)/(x1-x0)
end

function mapfrom01(x,y0,y1)
    return x*(y1-y0)+y0
end

function wavg(x0,x1,w=.5)
    return (1-w)*x0+w*x1
end


function step(x,x0=0)
    return x>x0 ? 1 : 0
end

function normalize(p)
    len = sqrt(p.x^2 + p.y^2)
    if len != 0
        p = Point(p.x / len, p.y / len)
    end
    return p
end

function rotate(p,angle)
    return Point(p.x*cos(angle)-p.y*sin(angle),p.x*sin(angle)+p.y*cos(angle))
end

function Base.:/(p1::Point, p2::Point)
    return Point(p1.x / p2.x, p1.y / p2.y)
end

function Base.:-(num::Number, p::Point)
    return Point(num - p.x, num - p.y)
end

function Base.:-(p1::Point, num::Number)
    return Point(p1.x - num, p1.y - num)
end

function Base.isless(num::Number, p::Point)
    return num < p.x && num < p.y
end

function Base.isless(p::Point, num::Number)
    return p.x < num && p.y < num
end

function Base.abs(p::Point)
    return Point(abs(p.x), abs(p.y))
end

function Base.isnan(p::Point)
    return isnan(p.x) || isnan(p.y)
end

function Base.min(p1::Point, p2::Point)
    return Point(min(p1.x, p2.x), min(p1.y, p2.y))
end

function Base.max(p1::Point, p2::Point)
    return Point(max(p1.x, p2.x), max(p1.y, p2.y))
end

function randF(min,max)
    return min+rand()*(max-min)
end

function randF(max)
    return randF(0,max)
end

function randI(min,max)
    return floor(randF(min,max))
end

function randI(max)
    return randI(0,max)
end

function linspace(start, stop, len)
    return range(start, stop=stop, length=len)
end

function meshgrid(x,y)
    return [Point(x,yi) for yi in y],[Point(xi,y) for xi in x]
end

function meshgrid(x_points::Int, y_points::Int, x_range::Tuple{Float64, Float64}, y_range::Tuple{Float64, Float64})
    x_lin = range(x_range[1], stop=x_range[2], length=x_points)
    y_lin = range(y_range[1], stop=y_range[2], length=y_points)
    
    x_grid = repeat(reshape(x_lin, 1, :), y_points, 1)
    y_grid = repeat(y_lin, 1, x_points)
    
    return x_grid, y_grid
end

function truncate(x)
    return x>0 ? floor(x) : ceil(x)
end

function Base.round(x,step)
    return round(x/step)*step
end

function parabola(x,root0=0,root1=1)
    a = 1 / ((0.5 * (root0 + root1) - root0) * (0.5 * (root0 + root1) - root1))
    return a*(x-root0)*(x-root1)/((root0-root1)^2)
end

function limit(x,x0,x1)
    return min(max(x,x0),x1)
end

function limitMag(x,x0)
    return limit(x,-x0,x0)
end

function cross(p1::Point,p2::Point)
    return p1.x*p2.y-p1.y*p2.x
end

function dot(p1::Point,p2::Point)
    return p1.x*p2.x+p1.y*p2.y
end

function sP(x,p) #Symmetric Power
    return abs(x)^p*sign(x)
end

function heading(p::Point)
    return atan(p.y,p.x)
end

function modulo(x,y=1)
    return x-y*floor(x/y)
end

function Base.round(x::Point; digits::Integer=0)
    return Point(round(x.x; digits=digits), round(x.y; digits=digits))
end


function smoothStep(x)
    if x<0
        return 0
    elseif x>1
        return 1
    else
        return 3*x^2-2*x^3
    end
end

# function scaledEase(x,c=0)
#     if c>0
#         return x/(c*(1-x)+1)
#     else
#         return 1-(1-x)/(-c*(x)+1)
#     end
# end

function scaledEase(x,c=0,x0=0,x1=1,y0=0,y1=1)
    x = mapto01(x,x0,x1)
    if c>0
        return mapfrom01(x/(c*(1-x)+1),y0,y1)
    else
        return mapfrom01(1-(1-x)/(-c*(x)+1),y0,y1)
    end
    # return mapfrom01(scaledEase(mapto01(x,x0,x1),c),y0,y1)
end

function scaledStep(x,c=0,x0=0,x1=1,y0=0,y1=1)
    x = mapto01(x,x0,x1)
    y =0
    if x>.5
        y= scaledEase(x,c,.5,1,.5,1)
    else
        y= scaledEase(x,c,.5,0,.5,0)
    end
    y = mapfrom01(y,y0,y1)
    # return mapfrom01(scaledStep(mapto01(x,x0,x1),c),y0,y1)
end

# function scaledStep(x,c=0)
#     if x>.5
#         # return mapfrom01(scaledEase(mapto01(x,.5,1),c),.5,1)
#         return scaledEase(x,c,.5,1,.5,1)
#     else
#         # return mapfrom01(scaledEase(mapto01(x,.5,0),c),.5,0)
#         return scaledEase(x,c,.5,0,.5,0)
#     end
# end

function limitL(p::Point,min,max)
    L=mag(p)
    L==0 && return p
    if L>max
        return p*max/L
    elseif L<min
        return p*min/L
    else
        return p
    end
end

function distance_to_segment(P, A, B)
    AB = B - A
    AP = P - A
    BP = P - B
    # Project point P onto line AB
    projection = A + dot(AP, AB) / dot(AB, AB) * AB
    # Check if projection lies within segment AB
    if dot(projection - A, projection - B) <= 0
        return mag(projection - P)
    else
        return min(mag(AP), mag(BP))
    end
end

function dist(p1::Point,p2::Point)
    return mag(p1-p2)
end

function Base.atan(p::Point)
    return atan(p.y,p.x)
end


horSymm(p) = Point(p.x, -p.y)
verSymm(p) = Point(-p.x, p.y)
diagSymm(p) = Point(p.y, p.x)
originSymm(p) = Point(-p.x, -p.y)
antiDiagSymm(p) = Point(-p.y, -p.x)


function angleBetweenVectors(v1, v2)
    dot_product = dot(v1, v2)
    mag_v1 = sqrt(dot(v1, v1))
    mag_v2 = sqrt(dot(v2, v2))
    
    cos_theta = dot_product / (mag_v1 * mag_v2)
    
    return acos(clamp(cos_theta, -1.0, 1.0)) # to handle possible floating-point inaccuracies
end

function signedAngleBetweenVectors(v1, v2)
    dot_product = dot(v1, v2)
    determinant = v1[1]*v2[2] - v1[2]*v2[1]
    
    mag_v1 = sqrt(dot(v1, v1))
    mag_v2 = sqrt(dot(v2, v2))
    
    cos_theta = dot_product / (mag_v1 * mag_v2)
    
    angle_radians = acos(clamp(cos_theta, -1.0, 1.0))
    
    # determine the sign based on the determinant
    return sign(determinant) * angle_radians
end

function vertexAngle(A::Point, B::Point, C::Point)
    BA = A-B
    BC = C-B

    angle_radians = signedAngleBetweenVectors(BA, BC)
    
    return angle_radians
end