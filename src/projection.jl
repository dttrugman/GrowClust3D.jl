### Helper Functions for Forward and Inverse Transformations ###

### Forward transform using transform object fproj ###
function lonlat2xypos(lons::Vector{Float64},lats::Vector{Float64},
    rotANG::Float64,fproj) # forward projection object

    # initialize
    n = length(lons)
    xx, yy = zeros(n), zeros(n)

    # project each lon lat point
    for ii in eachindex(xx)
        if rotANG == 0.0 # no rotation
            #xx[ii], yy[ii] = fproj([lons[ii] lats[ii]]) # v1.4 or lower
            xx[ii], yy[ii] = fproj(lons[ii], lats[ii]) # v1.5
        else # forward rotation after projection
            #xx0, yy0 = fproj([lons[ii] lats[ii]]) # v1.4 or lower
            xx0, yy0 = fproj(lons[ii], lats[ii]) # v1.5
            xx4 = xx0*cosd(rotANG) + yy0*sind(rotANG)
            yy4 = -xx0*sind(rotANG) + yy0*cosd(rotANG)
            xx[ii], yy[ii] = xx4, yy4
        end
    end

    # return
    return xx, yy

end

### Inverse transform using transform object iproj
function xypos2latlon(xx::Vector{Float64},yy::Vector{Float64},
    rotANG::Float64,iproj) # inverse projection object

    # initialize
    n = length(xx)
    lons, lats = zeros(n), zeros(n)

    # inverse project each lon lat point
    for ii in eachindex(xx)
        if rotANG == 0.0 # no rotation
            #lons[ii], lats[ii] = iproj([xx[ii] yy[ii]]) # v1.4 or lower
            lons[ii], lats[ii] = iproj(xx[ii], yy[ii]) # v1.5
        else # inverse rotation before inverse projection
            xx4, yy4 = xx[ii], yy[ii]
            xx0 = xx4*cosd(-rotANG) + yy4*sind(-rotANG) # negative to reverse
            yy0 = -xx4*sind(-rotANG) + yy4*cosd(-rotANG) # negative to reverse
            #lons[ii], lats[ii] = iproj([xx0 yy0]) # v1.4 or lower
            lons[ii], lats[ii] = iproj(xx0, yy0) # v1.5
        end
    end

    # return
    return lons, lats

end

### Function to setup map projection based on input parameters

function setup_projection(inpD,plon0,plat0,plat1,plat2)

    ## setup projection
    mapproj = inpD["proj"]
    rellipse = inpD["rellipse"]
    if mapproj in ["aeqd", "tmerc"]
        proj = Transformation("+proj=longlat +ellps=$rellipse",
            "+proj=$mapproj +ellps=$rellipse +lat_0=$plat0 +lon_0=$plon0 +units=km")
    elseif mapproj == "merc"
        proj = Transformation("+proj=longlat +ellps=$rellipse",
            "+proj=$mapproj +ellps=$rellipse +lat_ts=$plat0 +lon_0=$plon0 +units=km")
    elseif mapproj == "lcc"
        proj = Transformation("+proj=longlat +ellps=$rellipse",
        "+proj=$mapproj +ellps=$rellipse +lat_0=$plat0 +lon_0=$plon0 +lat_1= $plat1 +lat_2=$plat2 +units=km")
    else
        println("ERROR, map projection not defined! ", mapproj)
        exit()
    end

    # return projection object
    return proj

end

### Simple functions to compute distances between points on a map ####
#  > Inputs: lat1, lon1, lat2, lon2 (float or vector)
#  < Returns: distance in km  (float or vector)

# scalar to scalar
function map_distance(lat1::Float64,lon1::Float64,
    lat2::Float64,lon2::Float64,degkm=111.1949266)
lat0 = 0.5*(lat2+lat1)
return degkm*sqrt( (lat2-lat1)^2 + (cosd(lat0)*(lon2-lon1))^2 )
end

# scalar to vector
function map_distance(lat1::Float64,lon1::Float64,
    lat2::Vector{Float64},lon2::Vector{Float64},degkm=111.1949266)
lat0 = 0.5*(mean(lat2)+lat1)
return degkm*sqrt.( (lat2.-lat1).^2 .+ (cosd(lat0)*(lon2.-lon1)).^2 )
end

# vector to scalar
function map_distance(lat1::Vector{Float64},lon1::Vector{Float64},
    lat2::Float64,lon2::Float64,degkm=111.1949266)
lat0 = 0.5*(mean(lat1)+lat2)
return degkm*sqrt.( (lat2.-lat1).^2 .+ (cosd(lat0)*(lon2.-lon1)).^2 )
end

# vector to vector
function map_distance(lat1::Vector{Float64},lon1::Vector{Float64},
    lat2::Vector{Float64},lon2::Vector{Float64},degkm=111.1949266)
lat0 = 0.5*mean(lat2.+lat1)
return degkm*sqrt.( (lat2.-lat1).^2 .+ (cosd(lat0)*(lon2.-lon1)).^2 )
end

### Functions for computing cartestian distances

# vectors (x,y) x 2
function xydist(x1::Vector{Float64},y1::Vector{Float64},x2::Vector{Float64},y2::Vector{Float64})
    return sqrt.((x1.-x2).^2 .+ (y1.-y2).^2)
end

# vectors (x,y) with scalar reference
function xydist(x1::Vector{Float64},y1::Vector{Float64},x2::Float64,y2::Float64)
    return sqrt.((x1.-x2).^2 .+ (y1.-y2).^2)
end