### READ_GCINP reads the input file controlling script parameters
#
#  > Inputs: infile, relative path to file w/ algorithm parameters
#  < Returns: inpD, a dictionary with parsed control parameters
function read_gcinp(inpfile)
    ### open file
    counter = 0
    inpD = open(inpfile) do f
        
        # initialize dictionary, line counter
        inpD = Dict{String,Any}()

        # default table type
        inpD["ttabsrc"] = "trace"
        
        # loop over each line
        for line in eachline(f)
            
            # strip whitespace
            sline = strip(line)
            
            # check for empty or comment lines
            if length(sline) < 1
                continue
            elseif sline[1] == '*'
                continue
            else
                counter+=1
            end
            
            # parse initial lines
            if counter == 1
                inpD["evlist_fmt"] = parse(Int64,sline)
            elseif counter == 2
                inpD["fin_evlist"] = sline
            elseif counter == 3
                inpD["stlist_fmt"] = parse(Int64,sline)
            elseif counter == 4
                inpD["fin_stlist"] = sline
            elseif counter == 5
                data = split(sline)
                inpD["xcordat_fmt"] = parse(Int64,data[1])
                inpD["tdif_fmt"] = parse(Int64,data[2])
            elseif counter == 6
                inpD["fin_xcordat"] = sline
            elseif counter == 7
                inpD["ttabsrc"] = sline
            elseif counter == 8
                inpD["fin_vzmdl"] = sline
            elseif counter == 9
                inpD["fdir_ttab"] = sline
            elseif counter == 10
                data = split(sline)
                inpD["proj"] = data[1]
                inpD["rellipse"] = data[2]
                inpD["lon0"] = parse(Float64,data[3])
                inpD["lat0"] = parse(Float64,data[4])
                if inpD["proj"] == "lcc"
                    inpD["latp1"] = parse(Float64,data[5])
                    inpD["latp2"] = parse(Float64,data[6])
                    inpD["rotANG"] = parse(Float64,data[7])
                elseif length(data) >=5
                    inpD["rotANG"] = parse(Float64,data[5])
                else
                    inpD["rotANG"] = 0.0
                end
            elseif counter == 11
                data = split(sline)
                if inpD["ttabsrc"] == "trace"
                    inpD["vpvs_factor"] = parse(Float64,data[1])
                    inpD["rayparam_min"] = parse(Float64,data[2])
                else
                    inpD["vpvs_factor"] = parse(Float64,data[1])
                end
            elseif counter == 12
                data = split(sline)
                if inpD["ttabsrc"] == "trace"
                    inpD["tt_zmin"] = parse(Float64,data[1])
                    inpD["tt_zmax"] = parse(Float64,data[2])
                    inpD["tt_zstep"] = parse(Float64,data[3])
                elseif inpD["ttabsrc"] == "nllgrid"
                    inpD["tt_zmin"] = parse(Float64,data[1])
                    inpD["tt_zmax"] = parse(Float64,data[2])
                end
            elseif counter == 13
                data = split(sline)
                if inpD["ttabsrc"] == "trace"
                    inpD["tt_xmin"] = parse(Float64,data[1])
                    inpD["tt_xmax"] = parse(Float64,data[2])
                    inpD["tt_xstep"] = parse(Float64,data[3])
                    inpD["tt_ndim"] = 2
                elseif inpD["ttabsrc"] == "nllgrid"
                    if length(data) < 4
                        inpD["tt_xmin"] = parse(Float64,data[1])
                        inpD["tt_xmax"] = parse(Float64,data[2])
                        inpD["tt_ndim"] = 2
                    else
                        inpD["tt_xmin"] = parse(Float64,data[1])
                        inpD["tt_xmax"] = parse(Float64,data[2])
                        inpD["tt_ymin"] = parse(Float64,data[3])
                        inpD["tt_ymax"] = parse(Float64,data[4])
                        inpD["tt_ndim"] = 3
                    end
                end
            elseif counter == 14
                data = split(sline)
                inpD["rmin"] = parse(Float32,data[1])
                inpD["delmax"] = parse(Float64,data[2])
                inpD["rmsmax"] = parse(Float32,data[3])
            elseif counter == 15
                data = split(sline)
                inpD["rpsavgmin"] = parse(Float32,data[1])
                inpD["rmincut"] = parse(Float32,data[2])
                inpD["ngoodmin"] = parse(Int64,data[3])
                inpD["iponly"] = parse(Int64,data[4])
            elseif counter == 16
                data = split(sline)
                inpD["nboot"] = parse(Int64,data[1])
                inpD["nbranch_min"] = parse(Int64,data[2])
            elseif counter == 17
                inpD["fout_cat"] = sline
            elseif counter == 18
                inpD["fout_clust"] = sline
            elseif counter == 19
                inpD["fout_log"] = sline
            elseif counter == 20
                inpD["fout_boot"] = sline  
            else
                break
            end
            

        end # close loop over lines

        # check for too few lines
        if counter < 20
            println("Error, too few input lines.")
            println("Parsed inputs:\n",inpD)
            exit()
        else
            # return results
            return inpD
        end
    
    end # closes file
    
    # return
    return inpD

end

### INPUT_CHECK validates input file parameters
#
#  > Inputs: inpD, a dictionary with parsed control parameters
#  < Returns: input_ok, a boolean check
function check_gcinp(inpD)
    
    # inputs assumed ok until proven otherwise
    println("Checking input parameters...")
    input_ok = true

    # check vp/vs ratio
    if (inpD["vpvs_factor"] < 0.0)
        println("Input error, Vp/Vs vpvs_factor:")
        println("vpvs_factor") 
        input_ok = false
    end
    
    # check for travel time table source and related params
    if inpD["ttabsrc"] == "trace" # ray tracing

        # check min ray parameters plongcutP, plongcutS
        if  ((inpD["plongcutP"] < 0) | (inpD["plongcutS"] < 0))
            println("Input error, ray param cutoffs plongcutP, plongcut:")
            println(inpD["plongcutS"], ", ", inpD["plongcutS"])   
            input_ok = false
        end
    
    elseif inpD["ttabsrc"] == "nllgrid" # precomputed grids

        # check for valid directory storing these files
        if !(ispath(inpD["fdir_ttab"]))
            println("Input error, no path to travel time grids:")
            println(inpD["fdir_ttab"])
            input_ok = false
        end

    else # undefined
        println("Input error, undefined ttabsrc: ", inpD["ttabsrc"])
        input_ok=false
    end

    # check evlist and stlist format
    if !(inpD["evlist_fmt"] in [1,2,3])
        println( "Input error, unknown evlist format: ",
                inpD["evlist_fmt"])
        println("Must be either 1 or 2 or 3")
        println( "Please fix input file: ", inpfile)
        input_ok = false
     end
     if !(inpD["stlist_fmt"] in [1,2])
        println( "Input error, unknown stlist format: ",
                inpD["stlist_fmt"])
        println("Must be either 1 or 2")
        println( "Please fix input file: ", inpfile)
        input_ok = false
     end

    # check tdif format
    if ((inpD["tdif_fmt"] != 12) & (inpD["tdif_fmt"] != 21))
        println( "Input error, unknown tdif sign convention: ",
                inpD["tdif_fmt"])
        println( "   12: event 1 - event 2")
        println( "   21: event 2 - event 1")
        println( "Please fix input file: ", inpfile)
        input_ok = false
     end


    # check rmsmax and delmax
    if ((inpD["rmsmax"] <= 0.0) | (inpD["delmax"] <= 0.0) )
        println( "Input error, GrowClust params : rmsmax, delmax")
        println( "rmsmax, delmax")
        input_ok = false
    end

    # check iponly
    if ((inpD["iponly"] < 0) | (inpD["iponly"] > 2))
        println( "Input error, GrowClust params : iponly")
        println( "$iponly" )
        input_ok = false
    end

    # check nboot
    if ((inpD["nboot"] < 0) | (inpD["nboot"] > 10000)) 
         println( "Input error: nboot, maxboot")
         println( "$nboot, $maxboot")
         input_ok = false
    end
    
    # check map projections
    if !(inpD["proj"] in ["aeqd", "lcc", "merc", "tmerc"])
        println("parameter error: mapproj (not implemented)")
        println(inpD["proj"])
        input_ok = false
    end

    # check reference ellipses
    if !(inpD["rellipse"] in ["WGS84", "GRS80", "WGS72", "clrk80", "clrk66",
        "intl", "new_intl", "krass", "aust_SA", "airy", "bessel", "sphere"])
        println("parameter error: rellipse (not implemented)")
        println(inpD["rellipse"])
        input_ok = false
    end

    # return
    return input_ok
    
end

### CHECK_AUXPARAMS validates global parameters
#
#  > Inputs: global run parameters: hshiftmax, vshiftmax, rmedmax,
#       boxwid, nit, irelonorm, vzmodel_type
#  < Returns: params_ok, a boolean check
function check_auxparams(hshiftmax, vshiftmax, rmedmax,
    boxwid, nit, irelonorm, vzmodel_type)

# assume params ok unless problem is found
println("Checking auxiliary run parameters")
params_ok = true 

# check hshiftmax, vshiftmax
if ((hshiftmax <= 0.0) | (vshiftmax <= 0.0))  
    println("parameter error: hshiftmax, vshiftmax")
    println(hshiftmax)
    println(vshiftmax)
    params_ok = false
end

# check rmedmax
if (rmedmax < 0.0)  
    println("parameter error: rmedmax")
    println(rmedmax)
    params_ok = false
end 

# check boxwid, nit
if ((boxwid <= 0.0) | (nit >=200))  
    println("parameter error: boxwid, nit")
    println(boxwid, nit)
    params_ok = false
end

# check irelonorm
if ((irelonorm < 1) | (irelonorm > 3))  
    println("parameter error: irelonorm")
    println(irelonorm)
    params_ok = false
end

# check vz_model type (for 1D ray tracing)
if ((vzmodel_type < 1) | (vzmodel_type > 3))  
    println("parameter error: vzmodel_type")
    println(vzmodel_type)
    params_ok = false
end

# return
return params_ok                                                     


end

### Event Catalog Reader
#
#  > Inputs: name of event file, integer file format
#  < Returns: event dataframe
function read_evlist(evfile,evfmt)

    # define column names
    if evfmt == 1
        cols=["qyr","qmon","qdy","qhr","qmin","qsc",
            "qlat","qlon","qdep","qmag","eh","ez","rms","qid"]
    elseif evfmt == 2
        cols=["qyr","qmon","qdy","qhr","qmin","qsc",
            "qid","qlat","qlon","qdep","qmag"]
    elseif evfmt == 3
        cols=["qid","qyr","qmon","qdy","qhr","qmin","qsc",
            "qlat","qlon","qdep","rms","eh","ez","qmag"]
    end
    #ncols=length(cols)
    
    # read all columns, converting where possible
    df = DataFrame(readdlm(evfile,Float64),cols)
    
    # compute otime 
    df[!,:qsc] = Base.round.(df[!,:qsc],digits=3) # DateTimes have ms precision
    df[!,:qotime] = DateTime.(convert.(Int32,df[!,:qyr]),
        convert.(Int32,df[!,:qmon]),convert.(Int32,df[!,:qdy]),
        convert.(Int32,df[!,:qhr]),convert.(Int32,df[!,:qmin]),
        convert.(Int32,floor.(df[!,:qsc])), convert.(Int32,
            round.(1000.0*(df[!,:qsc].-floor.(df[!,:qsc]))))
        ) 

    # event ids and serial event number
    df[!,:qid] = convert.(Int64,df[!,:qid])
    df[!,:qix] = Vector{Int32}(1:nrow(df))
    
    # reformat for output
    outcols=["qix","qid","qotime","qmag","qlat","qlon","qdep"]
    select!(df,outcols)
    
    # return results
    return df
    
end

##### Simple Input Printing Function

function print_input(inpD)
    
    println("\nInput verified! Check results below:")
    println("====================================")
    @printf("Input files:\n")
    println("> Event list: ", inpD["fin_evlist"])
    println("> Event format: ", inpD["evlist_fmt"])
    println("> Station list: ", inpD["fin_stlist"])
    println("> Station format: ", inpD["stlist_fmt"])
    println("> Xcor dataset: ", inpD["fin_xcordat"])
    println("> Xcor format: ", inpD["xcordat_fmt"])
    println("> Tdif format: ", inpD["tdif_fmt"])
    println("> Travel time source: ", inpD["ttabsrc"])
    println("> Velocity/Grid model: ", inpD["fin_vzmdl"])
    println("Travel time grid directory: ",inpD["fdir_ttab"])
    @printf("Projection: %s %s %.6f %.6f %.2f\n",
        inpD["proj"],inpD["rellipse"],inpD["lon0"],inpD["lat0"],inpD["rotANG"])
    @printf("GrowClust parameters:\n")
    @printf("> rmin, delmax, rmsmax: %.3f %.1f %.3f\n",
        inpD["rmin"],inpD["delmax"],inpD["rmsmax"])
    @printf("> rpsavgmin, rmincut, ngoodmin, iponly: % .3f %.3f %d %d\n",
        inpD["rpsavgmin"],inpD["rmincut"],inpD["ngoodmin"],inpD["iponly"])
    @printf("> nboot, nbranch_min: %d %d\n",
        inpD["nboot"],inpD["nbranch_min"])
    @printf("Output files:\n")
    println("> Relocated catalog: ", inpD["fout_cat"])
    println("> Cluster file: ", inpD["fout_clust"])
    println("> Bootstrap file: ", inpD["fout_boot"])
    println("> Log file: ", inpD["fout_log"])

end

### Station List Reader
#
#  > Inputs: name of station file, integer file format
#  < Returns: event dataframe
function read_stlist(stfile,stfmt)

    # define column names
    if stfmt == 1
        cols=["sta", "slat", "slon"]
        dtypes = [String,Float64,Float64]
    elseif stfmt == 2
        cols=["sta", "slat", "slon", "selev"]
        dtypes = [String,Float64,Float64,Float64]
    end
    #ncols=length(cols)
    
    # read all columns
    df = DataFrame(readdlm(stfile,Any),cols)
    for (col, ctype) in zip(cols,dtypes)
        if ctype == String
            df[!,col] = string.(df[!,col]) 
        else
            df[!,col] = convert.(ctype,df[!,col])
        end
    end
    
    # convert elevations to km
    if "selev" in names(df)
        df[!,:selev] ./= 1000.0
    else # set to zero
        df[!,:selev] .= 0.0
    end
    
    # reformat for output, selecting unique rows
    outcols=["sta","slat","slon","selev"]
    select!(df,outcols)
    unique!(df) # removes duplicate rows
    
    # return results
    return df
    
end

### Xcorr Data Reader
#  > Inputs: input parameters, event dataframe, station dataframe
#  < Returns: xcor dataframe
#  note: original function for event / station data in latlon coords
function read_xcordata_lonlat(inpD,qdf,sdf)

    # unpack parameters
    xcfile = inpD["fin_xcordat"]
    xcfmt = inpD["xcordat_fmt"]
    tdiffmt = inpD["tdif_fmt"]
    rmincut = inpD["rmincut"]
    rpsavgmin = inpD["rpsavgmin"]
    rmingood = inpD["rmin"]
    ngoodmin =inpD["ngoodmin"]
    iponly = inpD["iponly"]
    delmax = inpD["delmax"]
    
    # only text file for now
    if xcfmt != 1
        println("XCOR FORMAT NOT IMPLEMENTED:",xcfmt)
        exit()
    end
    
    # read initial data
    cols = ["sta", "tdif", "rxcor", "iphase"]
    tmap = [String, Float32, Float32, Int8]
    xdf = DataFrame(readdlm(xcfile,Any,
        comments=true,comment_char='#', use_mmap=true),cols)
    Threads.@threads for ii in 1:length(cols)
        col = cols[ii]
        if ii == 1
            xdf[!,col] = string.(xdf[!,col])
        elseif ii == 4
            xdf[!,col] = ifelse.(xdf[!,col].=="P",Int8(1),Int8(2))
        else
            xdf[!,col] = convert.(tmap[ii],xdf[!,col])
        end
    end

   # add in event pairs
    xdf[!,:qid1] .= 0
    xdf[!,:qid2] .= 0    
    open(xcfile) do f
        ii = 0
        q1, q2 = 0, 0
        for line in eachline(f)
            if line[1]=='#'
                _, ev1, ev2, _ = split(line)
                q1 = parse(Int64,ev1)
                q2 = parse(Int64,ev2)
            else
                ii+=1
                xdf[ii,:qid1] = q1
                xdf[ii,:qid2] = q2
            end
        end
    end
            
    # Calculate Event Pair statistics
    transform!(groupby(xdf, [:qid1,:qid2]), :rxcor => mean => :gxcor)
    
    # Subset by quality
    if iponly==1
        xdf = xdf[(xdf[!,:rxcor].>=rmincut).&(
                   xdf[!,:gxcor].>=rpsavgmin).&(xdf[!,:iphase].==Int8(1)),:]
    else
        xdf = xdf[(xdf[!,:rxcor].>=rmincut).&(xdf[!,:gxcor].>=rpsavgmin),:]
    end
    
    # Merge Event Location
    xdf = innerjoin(xdf,qdf,on=:qid1=>:qid)
    DataFrames.rename!(xdf,:qlat=>:qlat1,:qlon=>:qlon1,:qix=>:qix1)
    xdf = innerjoin(xdf,qdf,on=:qid2=>:qid)
    DataFrames.rename!(xdf,:qlat=>:qlat2,:qlon=>:qlon2,:qix=>:qix2)
    
    # Centroid of event pair
    xdf[!,:qlat0] = 0.5*(xdf[!,:qlat1].+xdf[!,:qlat2])
    xdf[!,:qlon0] = 0.5*(xdf[!,:qlon1].+xdf[!,:qlon2])
    
    # Merge Station Location
    xdf = innerjoin(xdf,sdf,on=:sta)
    
    # Calculate distances
    xdf[!,:sdist] = map_distance(xdf[!,:qlat0],xdf[!,:qlon0],xdf[!,:slat],xdf[!,:slon])
    
    # Define "good" xcorr
    xdf[!,:igood] .= (xdf[!,:sdist].<=delmax).&(xdf[!,:rxcor].>=rmingood)
    transform!(groupby(xdf, [:qid1,:qid2]), :igood => sum => :ngood)
    
    # Keep only good pairs
    xdf = xdf[xdf[!,:ngood].>=ngoodmin,:]
    
    # redefine tdif to default tt2 - tt1
    if tdiffmt == 12
        xdf[!,:tdif].*=Float32(-1.0)
    end
    
    # Return       
    select!(xdf,["qix1","qid1","qix2","qid2","sta","tdif","rxcor",
            "slat","slon","sdist","igood","iphase"])
    return xdf
 
    
end

### Xcorr Data Reader
#  > Inputs: input parameters, event dataframe, station dataframe
#  < Returns: xcor dataframe
#  note: original function for event / station data in projected coords
function read_xcordata_proj(inpD,qdf,sdf)

    # unpack parameters
    xcfile = inpD["fin_xcordat"]
    xcfmt = inpD["xcordat_fmt"]
    tdiffmt = inpD["tdif_fmt"]
    rmincut = inpD["rmincut"]
    rpsavgmin = inpD["rpsavgmin"]
    rmingood = inpD["rmin"]
    ngoodmin =inpD["ngoodmin"]
    iponly = inpD["iponly"]
    delmax = inpD["delmax"]
    
    # only text file for now
    if xcfmt != 1
        println("XCOR FORMAT NOT IMPLEMENTED:",xcfmt)
        exit()
    end
    
    # read initial data
    cols = ["sta", "tdif", "rxcor", "iphase"]
    tmap = [String, Float32, Float32, Int8]
    xdf = DataFrame(readdlm(xcfile,Any,
        comments=true,comment_char='#', use_mmap=true),cols)
    Threads.@threads for ii in 1:length(cols)
        col = cols[ii]
        if ii == 1
            xdf[!,col] = string.(xdf[!,col])
        elseif ii == 4
            xdf[!,col] = ifelse.(xdf[!,col].=="P",Int8(1),Int8(2))
        else
            xdf[!,col] = convert.(tmap[ii],xdf[!,col])
        end
    end

   # add in event pairs
    xdf[!,:qid1] .= 0
    xdf[!,:qid2] .= 0    
    open(xcfile) do f
        ii = 0
        q1, q2 = 0, 0
        for line in eachline(f)
            if line[1]=='#'
                _, ev1, ev2, _ = split(line)
                q1 = parse(Int64,ev1)
                q2 = parse(Int64,ev2)
            else
                ii+=1
                xdf[ii,:qid1] = q1
                xdf[ii,:qid2] = q2
            end
        end
    end
            
    # Calculate Event Pair statistics
    transform!(groupby(xdf, [:qid1,:qid2]), :rxcor => mean => :gxcor)
    
    # Subset by quality
    if iponly==1
        xdf = xdf[(xdf[!,:rxcor].>=rmincut).&(
                   xdf[!,:gxcor].>=rpsavgmin).&(xdf[!,:iphase].==Int8(1)),:]
    else
        xdf = xdf[(xdf[!,:rxcor].>=rmincut).&(xdf[!,:gxcor].>=rpsavgmin),:]
    end
    
    # Merge Event Location
    xdf = innerjoin(xdf,qdf,on=:qid1=>:qid)
    DataFrames.rename!(xdf,:qX4=>:qX1,:qY4=>:qY1,:qix=>:qix1)
    xdf = innerjoin(xdf,qdf,on=:qid2=>:qid)
    DataFrames.rename!(xdf,:qX4=>:qX2,:qY4=>:qY2,:qix=>:qix2)
    
    # Centroid of event pair
    xdf[!,:qX4] = 0.5*(xdf[!,:qX1].+xdf[!,:qX2])
    xdf[!,:qY4] = 0.5*(xdf[!,:qY1].+xdf[!,:qY2])
    
    # Merge Station Location
    xdf = innerjoin(xdf,sdf,on=:sta)
    
    # Calculate distances
    #xdf[!,:sdist] = sqrt.((xdf.qX4.-xdf.sX4).^2 .+ (xdf.qY4.-xdf.sY4).^2)
    xdf[!,:sdist] = xydist(xdf.qX4, xdf.qY4, xdf.sX4, xdf.sY4)
    
    # Define "good" xcorr
    xdf[!,:igood] .= (xdf[!,:sdist].<=delmax).&(xdf[!,:rxcor].>=rmingood)
    transform!(groupby(xdf, [:qid1,:qid2]), :igood => sum => :ngood)
    
    # Keep only good pairs
    xdf = xdf[xdf[!,:ngood].>=ngoodmin,:]
    
    # redefine tdif to default tt2 - tt1
    if tdiffmt == 12
        xdf[!,:tdif].*=Float32(-1.0)
    end
    
    # Return       
    select!(xdf,["qix1","qid1","qix2","qid2","sta","tdif","rxcor",
            "sX4","sY4","sdist","igood","iphase"])
    return xdf
 
    
end