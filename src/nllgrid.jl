### Check projections in NonLinLoc Header Against Proj4
function check_proj(hdrD,inpD)
    
    # ok until proven otherwise
    proj_ok = true
    
    # check projection
    pmatches = [("AZIMUTHAL_EQUIDIST","aeqd",), # NLL and Proj4 matches
    ("LAMBERT","lcc"),("MERC","merc"), ("TRANS_MERC","tmerc")] 
    if !((hdrD["proj"], inpD["proj"]) in pmatches) 
        println("ERROR: PROJECTION TYPE MISMATCH")
        println(inpD["proj"], " & ", hdrD["proj"])
        proj_ok = false
    end

    # check reference ellipse
    rmatches = [
        ("WGS-84","WGS84"),("GRS-80","GRS80"),("WGS-72","WGS72"),
        ("Australian","aust_SA"),("Krasovsky","krass"),("International","new_intl"),
        ("Hayford-1909","intl"),("Clarke-1880","clrk80"),("Clarke-1866","clrk66"),
        ("Airy","airy"),("Bessel","bessel"),("Sphere","sphere"),
    ]
    if !((hdrD["ellipse"],inpD["rellipse"]) in rmatches)
        println("ERROR: PROJECTION ELLIPSOID MISMATCH")
        println(inpD["rellipse"], " & ", hdrD["ellipse"])
        proj_ok = false
    end

    # check lat0
    if inpD["lat0"] != hdrD["latORG"]
        println("ERROR: PROJECTION LAT0 MISMATCH")
        println(inpD["lat0"], " & ", hdrD["latORG"])
        proj_ok = false
    end
    
    # check lon0
    if inpD["lon0"] != hdrD["lonORG"]
        println("ERROR: PROJECTION LON0 MISMATCH")
        println(inpD["lon0"], " & ", hdrD["lonORG"])
        proj_ok = false
    end
    
    # check rotation angle
    if hdrD["rot_angle"] != inpD["rotANG"]
        println("ERROR: PROJECTION ROT_ANGLE MISMATCH")
        println(inpD["rotANG"], " & ", hdrD["rot_angle"])
        proj_ok = false    
    end
            
    # check parallels (LAMBERT ONLY)
    if hdrD["proj"] == "LAMBERT"
        if !((hdrD["parallel1"] == inpD["latp1"]) & (hdrD["parallel2"] == inpD["latp2"]))
            println("ERROR: PROJECTION PARALLEL MISMATCH")
            proj_ok = false
        end
    end

    # hope this stays at 0
    return proj_ok
end

### Read NonLinLoc Travel Time Table Header
#   input: path to header file (extension optional)
#   returns: dictionary of header information
function read_nll_head(hdrfile)
    
    # add extension
    if hdrfile[end-4:end] != ".hdr"
        hdrfile *= ".hdr"
    end
    
    # open file for reading
    hdrD = open(hdrfile) do f
        
        # initialize dictionary, line counter
        counter = 0
        hdrD = Dict{String,Any}()
        
        # loop over lines
        for line in eachline(f)
           
            # strip whitespace
            vals = split(strip(line))
            counter+=1
            
            # skip empty line
            if length(vals) < 1
                counter -=1
            
            # parse grid geometry
            elseif counter == 1
                hdrD["nX"] = parse(Int32,vals[1])
                hdrD["nY"] = parse(Int32,vals[2])
                hdrD["nZ"] = parse(Int32,vals[3])
                hdrD["xORG"] = parse(Float32,vals[4])
                hdrD["yORG"] = parse(Float32,vals[5])
                hdrD["zORG"] = parse(Float32,vals[6])
                hdrD["dX"] = parse(Float32,vals[7])
                hdrD["dY"] = parse(Float32,vals[8])
                hdrD["dZ"] = parse(Float32,vals[9])
                hdrD["gtype"] = vals[10]
                hdrD["dtype"] = vals[11]
                
            # parse transformation line if necessary
            elseif vals[1] in ["TRANS", "TRANSFORM"]
                hdrD["proj"] = vals[2]
                if hdrD["proj"] == "NONE"
                    continue
                elseif hdrD["proj"] in ["SDC", "SIMPLE"]
                    hdrD["latORG"] = parse(Float64,vals[4])
                    hdrD["lonORG"] = parse(Float64,vals[6])
                    hdrD["rot_angle"] = parse(Float64,vals[8])
                elseif hdrD["proj"] == "LAMBERT"
                    hdrD["ellipse"] = vals[4]
                    hdrD["latORG"] = parse(Float64,vals[6])
                    hdrD["lonORG"] = parse(Float64,vals[8])
                    hdrD["parallel1"] = parse(Float64,vals[10])
                    hdrD["parallel2"] = parse(Float64,vals[12])
                    hdrD["rot_angle"] = parse(Float64,vals[14])
                elseif hdrD["proj"] in ["TRANS_MERC", "AZIMUTHAL_EQUIDIST"]
                    hdrD["ellipse"] = vals[4]
                    hdrD["latORG"] = parse(Float64,vals[6])
                    hdrD["lonORG"] = parse(Float64,vals[8])
                    hdrD["rot_angle"] = parse(Float64,vals[10])
                end
             
            # station line
            else
                hdrD["station"] = vals[1]
                hdrD["staX"] = parse(Float32,vals[2])
                hdrD["staY"] = parse(Float32,vals[3])
                hdrD["staZ"] = parse(Float32,vals[4])
            end
            
            
        end
        
        # return parsed data
        return hdrD 
    end
    
    # return data as a dictionary
    return hdrD
end

###########################################

### Create Travel Time Interpolator
#   input: path to time grid file, associated header
#   returns: interpolation object
function make_nll_interp(grdfile,params)
    
    # extract grid parameters
    xmin, ymin, zmin = params["xORG"], params["yORG"], params["zORG"]
    dX, dY, dZ = params["dX"], params["dY"], params["dZ"]
    nX, nY, nZ = params["nX"], params["nY"], params["nZ"]
    xmax = xmin + dX*Float64(nX-1)
    ymax = ymin + dY*Float64(nY-1)
    zmax = zmin + dZ*Float64(nZ-1)

    # generate grid arrays
    xs = xmin:dX:xmax # not used in 2D
    ys = ymin:dY:ymax # y-->del in 2D
    zs = zmin:dZ:zmax # z-->dep in 2D
    
    # add extension
    if grdfile[end-4:end] != ".buf"
        grdfile *= ".buf"
    end
    
    # read binary data (Float32)
    ndata = nX*nY*nZ
    data = Array{Float32}(undef, ndata)
    read!(grdfile,data)

    # apply Vp/Vs scaling if necessary (e.g., convert P grid to S)
    if params["vscale"] != 1.0
        data *= Float32(params["vscale"])
    end

    # reshape into multidimensional grid: either (R,Z) or (X,Y,Z)
    if params["gtype"] == "TIME2D"
        data = permutedims(reshape(data,nZ,nY),(2,1)) # output is R/Z order, Float32
        
    else
        data = permutedims(reshape(data,nZ,nY,nX),(3,2,1)) # output is X/Y/Z, Float32
    end
    
    # subset grid boundaries (can save memory)
    if "xbounds" in keys(params)
        if params["gtype"] == "TIME2D"
            ymin, ymax = params["xbounds"]
            idxY = findfirst(x->x>=ymin,ys)
            ymin = ys[idxY] # first on-grid point
            ymax = min(ymax,ys[end]) # last on-grid point
            zmin, zmax = params["zbounds"]
            idxZ = findfirst(x->x>=zmin,zs)
            zmin = zs[idxZ] # first on-grid point
            zmax = min(zmax,zs[end]) # last on-grid point
            ys = ymin:dY:ymax # updates grid
            zs = zmin:dZ:zmax # updates grid
            jdxY = idxY+length(ys)-1
            jdxZ = idxZ+length(zs)-1
            data = data[idxY:jdxY,idxZ:jdxZ] # slice data matrix
        else # X, Y, Z
            xmin, xmax = params["xbounds"]
            idxX = findfirst(x->x>=xmin,xs)
            xmin = xs[idxX] # first on-grid point
            xmax = min(xmax,xs[end]) # last on-grid point
            ymin, ymax = params["ybounds"]
            idxY = findfirst(x->x>=ymin,ys)
            ymin = ys[idxY] # first on-grid point
            ymax = min(ymax,ys[end]) # last on-grid point
            zmin, zmax = params["zbounds"]
            idxZ = findfirst(x->x>=zmin,zs)
            zmin = zs[idxZ] # first on-grid point
            zmax = min(zmax,zs[end]) # last on-grid point
            xs = xmin:dX:xmax # updates grid
            ys = ymin:dY:ymax # updates grid
            zs = zmin:dZ:zmax # updates grid
            jdxX = idxX+length(xs)-1
            jdxY = idxY+length(ys)-1
            jdxZ = idxZ+length(zs)-1
            data = data[idxX:jdxX,idxY:jdxY,idxZ:jdxZ] # slice data matrix
        end
    end

    # extrapolation options for shallow depths
    if !("shallowmode" in keys(params))
        params["shallowmode"] = "throw"
    end
    if params["shallowmode"] == "flat"
        exconZ = Flat()
    elseif params["shallowmode"] == "reflect"
        exconZ = Reflect()
    elseif params["shallowmode"] == "linear"
        exconZ = Linear()
    else
        exconZ = Throw()
    end
    
    
    # define and interpolation object
    if !("interp_mode" in keys(params))
        params["interp_mode"] = "linear" # default to this
    end
    if params["gtype"] == "TIME2D" # 2D grid
        if params["interp_mode"] == "linear"
            return LinearInterpolation( (ys,zs), data, 
            extrapolation_bc=( (Line(),Line()), (exconZ,Line()) ) )
        elseif params["interp_mode"] == "cubic"
            return CubicSplineInterpolation((ys,zs), data,
            extrapolation_bc=( (Line(),Line()), (exconZ,Line()) ) )
        else
            println("Error, undefined interpolation: ",params["interp_mode"])
            return Nothing
        end
    else # 3D grid
        if params["interp_mode"] == "linear"
            return LinearInterpolation( (xs,ys,zs), data,
            extrapolation_bc=( (Line(),Line()), (Line(),Line()), (exconZ,Line()) ) )
        elseif params["interp_mode"] == "cubic"
            return CubicSplineInterpolation((xs,ys,zs), data,
            extrapolation_bc=( (Line(),Line()), (Line(),Line()), (exconZ,Line()) ) )
        else
            println("Error, undefined interpolation: ",params["interp_mode"])
            return Nothing
        end 
    end
end


