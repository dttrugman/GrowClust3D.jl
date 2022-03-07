### Script for Running GrowClust: Julia

### Import packages

# External Packages
using Printf
using DataFrames
using Random
using Dates
using Proj4: Transformation, inv

# GrowClust3D
using GrowClust3D

### Define Algorithm Parameters

# ------- GrowClust algorithm parameters -------------------------------
const distmax = 5.0          # maximum catalog(input) distance to join clusters (km)
const distmax2 = 3.0         # maximum relocated distance to join clusters (km)
const hshiftmax = 2.0        # maximum permitted horizontal cluster shifts (km)
const vshiftmax = 2.0        # maximum permitted vertical cluster shifts (km)
const rmedmax = Float32(0.05)         # maximum median absolute tdif residual to join clusters
const maxlink = 10           # use 10 best event pairs to relocate (optimize later...)
const nupdate = 1000         # update progress every nupdate pairs - NEW
   
# ------- Relative Relocation subroutine parameters -------------
const boxwid = 3. # initial "shrinking-box" width (km)
const nit = 15 # number of iterations
const irelonorm = 1 # relocation norm (L1 norm=1, L2 norm=2, 3=robust L2)
const tdifmax = Float32(30.) # maximum differential time value allowed (for error-checking on xcor data input)
const torgdifmax = Float32(10.0) # maximum origin time adjustment (for robustness)
   
# -------- Bootstrap resampling parameters -------------------
const iseed = 0 # random number seed
#  1: Resample each event pair independently (each pair always has the same # of picks in each resample)'
#  2: Resample the entire data vectors at once (event pairs may have different # of picks in each resample)
        
# ------- Velocity model parameters  -------------------------
const vzmodel_type = 1 # velocity model type: 1 = flat earth, (Z,Vp,Vs)
                 #                  or  2 = radial, (R,Vp,Vs): 
                 #                note: option 2 has not been extensively tested
const shallowmode = "flat" # option for how to treat shallow seismicity for 1D ray tracing
                       # flat treats negative depths as zero depth
                       # throw makes an error for negative depths
                       # linear does linear interpolation to negative depths
                       # reflect treats negative depths as equivalent to -depth

# ------- Geodetic parameters -------------
const degkm = 111.1949266 # for simple degree / km conversion
const erad = 6371.0 # average earth radius in km, WGS-84

### Read Input File

# read input file
println("\nReading input file: ")
infile_ctl = ARGS[1]
inpD = read_gcinp(infile_ctl)

### Update fields for ray tracing

# set default Vp/Vs
if inpD["vpvs_factor"] < 0.01
    inpD["vpvs_factor"] = sqrt(3.0)
end

# only needed for ray tracing
if inpD["ttabsrc"] == "trace" # only needed for ray tracing

    # define minimum ray parameter
    if inpD["rayparam_min"] < 0.0
        inpD["plongcutP"] = 1.0/7.8 # assumes a sub-moho Vp of 7.8 km/s
        inpD["plongcutS"] = sqrt(3)/7.8 # Vs = Vp/sqrt(3) by default
    else
        inpD["plongcutP"] = inpD["rayparam_min"] # user defined value
        inpD["plongcutS"] = inpD["rayparam_min"]*inpD["vpvs_factor"] # user defined value
    end

    # interpolation spacing
    inpD["itp_dz"] = inpD["tt_zstep"]
end


### Check Input parameters
input_ok = check_gcinp(inpD)
if input_ok
    println("Input parameters are ok!")
else
    println("ERROR: FIX INPUT PARAMETERS")
    exit()
end

# make output paths
println("Assigning output directories:")
for fkey in collect(keys(inpD))
    if startswith(fkey,"fout") # output catalogs, etc
        if inpD[fkey] != "NONE"
            outdir = join(split(inpD[fkey],"/")[1:end-1],"/")
            println(fkey, " => ",outdir)
            mkpath(outdir)
        end
    end
end
if (inpD["ttabsrc"] == "trace")&(inpD["fdir_ttab"]!="NONE")
    mkpath(inpD["fdir_ttab"]) # to store output tables
end


### Check Auxiliary Parameters
params_ok = check_auxparams(hshiftmax, vshiftmax, rmedmax,
        boxwid, nit, irelonorm, vzmodel_type)
if input_ok
    println("Auxiliary parameters are ok!")
else
    println("ERROR: FIX AUXILIARY PARAMETERS")
    exit()
end

### Print Input Parameters
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

### Read Catalog
println("\nReading event list:")
@time qdf = read_evlist(inpD["fin_evlist"],inpD["evlist_fmt"])
qid2qnum = Dict(zip(qdf.qid,qdf.qix))# maps event id to serial number
mlat0, mlon0 = median(qdf.qlat), median(qdf.qlon)
mlat01, mlat99 = percentile(qdf.qlat,1.0), percentile(qdf.qlat,99.0)
mlon01, mlon99 = percentile(qdf.qlon,1.0), percentile(qdf.qlon,99.0)
println("Median event location: $mlon0 $mlat0")
println("Longitude Percentiles (1%, 99%): $mlon01 $mlon99")
println("Latitude Percentiles (1%, 99%): $mlat01 $mlat99")

# validate projection
plat0, plon0 = inpD["lat0"], inpD["lon0"]
println("Projection origin: $plon0 $plat0")
if inpD["proj"] == "lcc"
    plat1, plat2 = inpD["latp1"], inpD["latp2"]
    println("Lambert Projection Parallels  : $plat1 $plat2")
end
if !(mlon01 < plon0 < mlon99)
    println("PROJECTION ORIGIN NOT ALIGNED WITH SEISMICITY:")
    exit()
elseif !(mlat01 < plat0 < mlat99)
    println("PROJECTION ORIGIN NOT ALIGNED WITH SEISMICITY:")
    exit()
end

### Read Stations
print("\nReading station list")
@time sdf = read_stlist(inpD["fin_stlist"],inpD["stlist_fmt"])
min_selev = minimum(sdf.selev)
max_selev = maximum(sdf.selev)
mean_selev = mean(sdf.selev)
@printf("station elevation (min,mean,max): %.3fkm %.3fkm %.3fkm\n",
    min_selev, mean_selev, max_selev)
sta2elev = Dict(zip(sdf.sta,sdf.selev)) # map station to elevation

#### Map Projection

## setup projection
const mapproj = inpD["proj"]
const rellipse = inpD["rellipse"]
const rotANG = inpD["rotANG"]
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
iproj = inv(proj) # inverse transformation

# project station coordinates, including rotation
sX4, sY4 = lonlat2xypos(sdf.slon,sdf.slat,rotANG,proj)
sdf[!,:sX4] .= sX4
sdf[!,:sY4] .= sY4

# project event coordinates
qX4, qY4 = lonlat2xypos(qdf.qlon,qdf.qlat,rotANG,proj)
qdf[!,:qX4] .= qX4
qdf[!,:qY4] .= qY4
println("\nProjected event list:")
show(qdf)
println()

### Read Xcor Data
println("\nReading xcor data") # in projected coordinates
@time xdf = read_xcordata_proj(inpD,qdf[!,[:qix,:qid,:qX4,:qY4]],sdf[!,[:sta,:sX4,:sY4]])
show(xdf)
println()
nbad = sum(abs.(xdf.tdif).>tdifmax)
if nbad > 0 # validate differential times
    println("Error: bad input differential times, some larger than tdifmax=$tdifmax")
    println("Fix input file or adjust tdifmax parameter.")
    ibad = abs.(xdf.tdif).>tdifmax
    show(xdf[ibad,:])
    exit()
end

### Gather Unique Station/Phase combos
usta = unique(xdf[!,"sta"])
nstaU = length(usta)
println("$nstaU unique stations used:\n",usta,"\n")
xdf[!,:itab] .= Int16(0)
staIDX = Dict(zip(usta,1:length(usta)))
for ii = 1:nrow(xdf) # table index: first nstaU are for P-waves, next are for S-waves
    xdf[ii,:itab] = staIDX[xdf[ii,:sta]] + nstaU*(xdf[ii,:iphase]-1)
end
ntab = 2*nstaU
println("Updated with table index:")
show(xdf[!,[:qix1,:qix2,:sta,:tdif,:rxcor,:iphase,:itab]])
println()

### Subset station list to get min/max values
subset!(sdf, :sta => ByRow(x -> x in usta))
minSX, maxSX = minimum(sdf.sX4), maximum(sdf.sX4)
minSY, maxSY = minimum(sdf.sY4), maximum(sdf.sY4)
minSR, maxSR = minimum(xdf.sdist), maximum(xdf.sdist)
@printf("min and max staX: %.4f %.4f\n", minSX, maxSX)
@printf("min and max staY: %.4f %.4f\n", minSY, maxSY)
@printf("min and max staR: %.4f %.4f\n", minSR, maxSR)

###### Ray-Tracing Section ##########


######## Pre-computed NonLinLoc Grids ########
if inpD["ttabsrc"] == "nllgrid"

    # store all tables here
    ttLIST = [] # stations for P, then stations for S
    phases = ["P","S"] # make P and S tables

    # loop over phases
    println("\nREADING TRAVEL TIME GRIDS:")
    total_time = @elapsed for (iphase, phase) in enumerate(phases)

        # loop over stations 
        for ista = 1:nstaU

            # filename
            sta = usta[ista] # model.phase.sta.time
            fname = inpD["fin_vzmdl"] * "." * phase * "." * sta *".time"
            println("Working on: $fname")

            # get header
            grdparams = read_nll_head(inpD["fdir_ttab"]*fname)

            # check projection
            proj_ok = check_proj(grdparams, inpD)
            if !(proj_ok)
                println("ERROR: PROJECTION MISMATCH: $fname")
                exit()
            end

            # check bounds
            if grdparams["gtype"] == "TIME2D"
                gmaxR = grdparams["yORG"] + grdparams["dY"]*Float64(grdparams["nY"]-1)
                if gmaxR < maxSR
                    println("ERROR: maximum station distance too large for grid!")
                    println("max station distance: $maxSR, max grid distance: $gmaxR")
                    exit()
                end
            else
                gminX = grdparams["xORG"]
                gmaxX = grdparams["xORG"] + grdparams["dX"]*Float64(grdparams["nX"]-1)
                gminY = grdparams["yORG"]
                gmaxY = grdparams["yORG"] + grdparams["dY"]*Float64(grdparams["nY"]-1)
                if ((gminX > minSX) | (gmaxX < maxSX) | (gminY > minSY) | (gmaxY < maxSY))
                    println("ERROR: Stations outside travel time grid!")
                    @printf("min and max staX: %.4f %.4f\n", minSX, maxSY)
                    @printf("min and max staY: %.4f %.4f\n", minSY, maxSY)
                    @printf("min and max gridX: %.4f %.4f\n",gminX,gmaxX)
                    @printf("min and max gridY: %.4f %.4f\n",gminY,gmaxY)
                    exit()
                end
            end
            
            # update header
            grdparams["interp_mode"] = "linear"
            grdparams["xbounds"] = [inpD["tt_xmin"],inpD["tt_xmax"]]
            if "tt_ymax" in keys(inpD)
                grdparams["ybounds"] = [inpD["tt_ymin"], inpD["tt_ymax"]]
            else
                grdparams["ybounds"] = [inpD["tt_xmin"], inpD["tt_xmax"]]
            end
            grdparams["zbounds"] = [inpD["tt_zmin"],inpD["tt_zmax"]]
            grdparams["shallowmode"] = shallowmode

            # get interpolator
            itpL = make_nll_interp(inpD["fdir_ttab"]*fname,grdparams)
            
            # update
            push!(ttLIST,itpL)

        end # end loop over stations
        
    end

    # finalize travel time tables
    const ttTABs = [deepcopy(xx) for xx in ttLIST] # static version
    const pTT = ttTABs[1] # p-wave example
    const sTT = ttTABs[Int64(ntab/2)+1] # s-wave example

else # Placeholder for now
    println("Travel time calculation mode not yet implemented: ", inpD["ttabsrc"])
    println("Ending program.")
    exit()
end

##### Test Interpolation Routines ####

# define test distances, depth
println("\nTesting travel time interpolation (P and S waves):")
test_dists = collect(range(0.0,60.0,step=5.0))
npts = length(test_dists)
test_depth = 10.0
     
# calculate travel times
test_ttP = pTT.(test_dists, test_dists, test_depth)
test_ttS = sTT.(test_dists, test_dists, test_depth)

# print results
println(typeof(test_ttP[1]), " ", typeof(test_ttS[1]))
for ii in 1:npts     
@printf("XY position: %4.0f %4.0f, Depth: %3.0f",
    test_dists[ii],test_dists[ii],test_depth)
@printf(" --> P and S travel times: %5.2f, %5.2f\n",
    test_ttP[ii],test_ttS[ii])
end

#### Validate Event Depths and Travel Time Tables

println("\nChecking event depths...")

# print event depths
const min_qdep = minimum(qdf.qdep)
const max_qdep = maximum(qdf.qdep)
@printf("min and max event depth: %.3fkm %.3fkm\n",min_qdep,max_qdep)

# print table depths
@printf("min and max table depth: %.3fkm %.3fkm\n",inpD["tt_zmin"],inpD["tt_zmax"])

# implement warnings and checks
if (min_qdep < inpD["tt_zmin"])
    println("WARNING: min event depth < min table depth")
end
if (inpD["ttabsrc"]=="trace")
    if (min_qdep < z_s0[1]) # only checked for ray-trace
        println("WARNING: min event depth < min vzmodel depth")
        println("Results may be inaccurate near free-surface...")
        #exit() # allow this, but warn user (should be ok if depth is near 0)
    end
end
if (max_qdep > inpD["tt_zmax"]) # note tt_zmax is >= vzmax
    println("ERROR: max event depth > max table / velocity model depth")
    exit()
end

println("Done.")

#########################################################################

############# Main Clustering Loop: Including Bootstrapping ##############

# define event-based output arrays
const nq = Int32(nrow(qdf))
revids = qdf[:,:qid]
rlats = qdf[:,:qlat]
rlons = qdf[:,:qlon]
rXs = qdf[:,:qX4]
rYs = qdf[:,:qY4]
rdeps = qdf[:,:qdep]
rorgs = zeros(Float32,nq) # origin time adjust
rcids = Vector{Int32}(1:nq) # initialize each event into one cluster
npair = -1

# Setup bootstrapping matrices
if inpD["nboot"] > 0
    blatM = repeat(qdf.qlat,1,inpD["nboot"])
    blonM = repeat(qdf.qlon,1,inpD["nboot"])
    bdepM = repeat(qdf.qdep,1,inpD["nboot"])
    borgM = zeros(Float32,(nq,inpD["nboot"]))
    bnbM = repeat(Vector{Int32}(1:nq),1,inpD["nboot"])
end

# base xcor dataframe to sample from
xdf00 = select(xdf,[:qix1,:qix2,:tdif,:itab,:igood])
xdf00[!,:gxcor] = ifelse.(xdf.igood.>0,xdf.rxcor,Float32(0.0)) # xcor with bad values zeroed
#show(xdf00)

# sampling vector
if nrow(xdf00)<typemax(Int32)
    const nxc = Int32(nrow(xdf00))
    const ixc = Vector{Int32}(1:nxc)
else
    const nxc = nrow(xdf00)
    const ixc = Vector{Int64}(1:nxc)
end


#### loop over each bootstrapping iteration
println("\n\n\nStarting relocation estimates.")
@time for ib in 0:inpD["nboot"]    
    
    # timer for this thread
    wc = @elapsed begin
    
    # bootstrapping: resample data before run
    Random.seed!(iseed + ib) # different for each run
    wc2 = @elapsed begin
    println("Initializing xcorr data and event pairs. Bootstrap iteration: $ib")
    if ib > 0 # sample with replacement from original xcorr array
        isamp = sort(sample(ixc,nxc,replace=true)) # sorted to keep evpairs together
        rxdf = xdf00[isamp,:]
    else     # do not resample on first run
        rxdf = xdf00
    end

    # assign xcorr index array
    if nxc < typemax(Int32)
        rxdf[!,:ixx] = Vector{Int32}(1:nxc)
    else
        rxdf[!,:ixx] = Vector{Int64}(1:nxc)
    end

    # compile event pair arrays
    bpdf = combine(groupby(rxdf[!,Not([:itab,:tdif])],[:qix1,:qix2]),
        :gxcor=>sum=>:rfactor,:ixx=>first=>:ix1,:ixx=>last=>:ix2,:igood=>sum=>:ngood)
    
    # sort pairs (note, resampled pairs may not have ngoodmin tdifs)
    if ib > 0
        bpdf = bpdf[bpdf.ngood.>=inpD["ngoodmin"]-2,[:qix1,:qix2,:rfactor,:ix1,:ix2]] # for robustness
    else
        select!(bpdf,[:qix1,:qix2,:rfactor,:ix1,:ix2])
    end
    sort!(bpdf,:rfactor,rev=true) # so best pairs first
    #show(bpdf)
    end # ends elapsed time for setup
    println("\nDone, elapsed time = $wc2")
    
    # run clustering
    brXs, brYs, brdeps, brorgs, brcids, bnb = clustertree_3D(
        bpdf.qix1, bpdf.qix2, bpdf.ix1, bpdf.ix2, 
        rxdf.tdif, rxdf.itab,
        qdf.qX4, qdf.qY4, qdf.qdep,
        ttTABs, nit, boxwid, irelonorm,
        inpD["rmsmax"],rmedmax,distmax,distmax2,
        hshiftmax,vshiftmax,torgdifmax,nupdate,maxlink)

    # inverse projection back to map coordinates
    brlons, brlats = xypos2latlon(brXs,brYs,rotANG,iproj)
    
    # save output
    if ib > 0
        blatM[:,ib] .= brlats
        blonM[:,ib] .= brlons
        bdepM[:,ib] .= brdeps
        borgM[:,ib] .= brorgs
        bnbM[:,ib] .= bnb
    else
        rXs .= brXs
        rYs .= brYs
        rlats .= brlats
        rlons .= brlons
        rdeps .= brdeps
        rorgs .= brorgs
        rcids .= brcids
        global npair = nrow(bpdf)
    end
        
    # completion
    end # ends the wall clock

end

################################################################

### Finalize Clustering Trees ###

println("\nFinalizing clustering trees...")

# temporary dataframe with event and cluster number
tdf = DataFrame("enum"=>qdf.qix,"cnum"=>rcids)

# compute nbranch
transform!(groupby(tdf, :cnum), nrow => :nb)

# unique clusters, sort by nbranch
select!(tdf,[:cnum,:nb])
unique!(tdf)
sort!(tdf,[:nb,:cnum],rev=[true,false])

# assign new cluster ids, starting with largest
tdf[!,:cid] =range(1,nrow(tdf),step=1)

# add back in event ids, reorder by event index
tdf = innerjoin(DataFrame("enum"=>qdf.qix,"cnum"=>rcids),tdf,on=:cnum)
sort!(tdf,:enum)
nclust = maximum(tdf.cid)
#show(tdf) 

# update cluster ids for all events
rcids00 = copy(rcids)
rnbranch = tdf.nb
rcids = tdf.cid

# finalize t-arrays
tlons, tlats = zeros(nclust), zeros(nclust) 
tdeps, torgs = zeros(nclust), zeros(Float32,nclust)
tnbranch = zeros(Int64,nclust)
for icc in 1:nclust
    idx = findall(rcids.==icc)
    tlons[icc] = mean(rlons[idx])
    tlats[icc] = mean(rlats[idx])
    tdeps[icc] = mean(rdeps[idx])
    torgs[icc] = mean(rorgs[idx])
    tnbranch[icc] = length(idx)
end

# compute cluster counts
ntree2 = sum(tnbranch.>=2)
ntree5 = sum(tnbranch.>=5)
ntree10 = sum(tnbranch.>=10)
ntree20 = sum(tnbranch.>=20)
ntree50 = sum(tnbranch.>=50)
ntree100 = sum(tnbranch.>=100)

### Relocated DataFrame

# initialize
println("\nFinalizing relocated dataset:")
rdf = DataFrame("enum"=>qdf.qix,"evid"=>revids,
    "rlat"=>rlats,"rlon"=>rlons,"rdep"=>rdeps,"rtim"=>rorgs,
    "rcid"=>rcids,"rnb"=>rnbranch,"rX"=>rXs,"rY"=>rYs)

# compute relocated
nreloc = sum(rdf.rnb .> 1)
println("Relocated: ",nreloc)

# adjusted otime
rdf[!,:rot] = qdf[!,:qotime] .+ Nanosecond.(round.(rorgs*1.0e9))
show(rdf)
println()


### Compute Misfits - w/ otime adjustment ###

println("\nComputing misfits...")

# select columns
resdf = select(xdf,[:qid1,:qid2,:tdif,:itab,:rxcor,:sdist])

# merge with event data, renaming columns to specify 1/2
resdf = innerjoin(resdf,rdf[!,[:enum,:evid,:rX,:rY,:rdep,:rtim]],on=:qid1=>:evid)
DataFrames.rename!(resdf,:enum=>:qnum1,:rX=>:qX1,:rY=>:qY1,:rdep=>:qZ1,:rtim=>:qtim1)
resdf = innerjoin(resdf,rdf[!,[:enum,:evid,:rX,:rY,:rdep,:rtim]],on=:qid2=>:evid)
DataFrames.rename!(resdf,:enum=>:qnum2,:rX=>:qX2,:rY=>:qY2,:rdep=>:qZ2,:rtim=>:qtim2)

# here, compute misfits for relocated event pairs and good differential times (for consistency with f90 codes)
resdf = resdf[(rcids[resdf.qnum1] .== rcids[resdf.qnum2]) .& ( # only event pairs in same cluster
    resdf[!,:sdist].<=inpD["delmax"]).&(resdf[!,:rxcor].>=inpD["rmin"]),:]

# compute predicted travel times
resdf[!,:pdif] .= 0.0
for ii = 1:nrow(resdf)
    resdf[ii,:pdif] = ttTABs[resdf[ii,:itab]](resdf[ii,:qX2],resdf[ii,:qY2],resdf[ii,:qZ2]) -
        ttTABs[resdf[ii,:itab]](resdf[ii,:qX1],resdf[ii,:qY1],resdf[ii,:qZ1]) +
        resdf[ii,:qtim2] - resdf[ii,:qtim1]
end

# P vs S
ipp = (resdf[!,:itab].<=ntab/2)
npp = sum(ipp)
iss = .!ipp
nss = sum(iss)

# RMS
rmsP = evalrms(resdf[ipp,:tdif].-resdf[ipp,:pdif])
rmsS = evalrms(resdf[iss,:tdif].-resdf[iss,:pdif])

# Mean signed residuals
msresP = mean(resdf[ipp,:tdif].-resdf[ipp,:pdif])
msresS = mean(resdf[iss,:tdif].-resdf[iss,:pdif])

# print
println("P-wave RMS: $rmsP")
println("Mean signed residual: $msresP")
println("Phases used: $npp")
println("S-wave RMS: $rmsS")
println("Mean signed residual: $msresS")
println("Phases used: $nss")

# show
select!(resdf,[:qnum1,:qnum2,:qid1,:qid2,:tdif,:itab,:pdif])

### Compute Event-based Stats

# arrays to store stats for each event
qnpair = zeros(Int64,nq)
qndiffP, qndiffS = zeros(Int64,nq), zeros(Int64,nq)
qsseP, qsseS = zeros(nq), zeros(nq)

# event 1 in pair, P-wave
for subdf in groupby(resdf[ipp,:],:qnum1)
    qix=subdf[1,:qnum1]
    qndiffP[qix] += nrow(subdf)
    qsseP[qix] += sum((subdf[!,:tdif].-subdf[!,:pdif]).^2)
end

# event 1 in pair, S-wave
for subdf in groupby(resdf[iss,:],:qnum1)
    qix=subdf[1,:qnum1]
    qndiffS[qix] += nrow(subdf)
    qsseS[qix] += sum((subdf[!,:tdif].-subdf[!,:pdif]).^2)
end

# event 2 in pair, P-wave
for subdf in groupby(resdf[ipp,:],:qnum2)
    qix=subdf[1,:qnum2]
    qndiffP[qix] += nrow(subdf)
    qsseP[qix] += sum((subdf[!,:tdif].-subdf[!,:pdif]).^2)
end

# event 2 in pair, S-wave
for subdf in groupby(resdf[iss,:],:qnum2)
    qix=subdf[1,:qnum2]
    qndiffS[qix] += nrow(subdf)
    qsseS[qix] += sum((subdf[!,:tdif].-subdf[!,:pdif]).^2)
end

# event pairs
for subdf in groupby(resdf,[:qnum1,:qnum2])
    qix1, qix2 = values(subdf[1,[:qnum1,:qnum2]])
    qnpair[qix1] += 1
    qnpair[qix2] += 1
end

# combined data --> RMS values
qrmsP = ifelse.(qndiffP.>0,sqrt.(qsseP./qndiffP),NaN64)
qrmsS = ifelse.(qndiffS.>0,sqrt.(qsseS./qndiffS),NaN64)

### Write Log / Statistics
println("Done.\n\nRUN SUMMARY TO FOLLOW:\n")

# Run Parameters
@printf("********************** Input Parameters ***********************\n")
@printf("     control file:   %s\n", infile_ctl)
@printf("       event list:   %s\n", inpD["fin_evlist"])
@printf("     station list:   %s\n", inpD["fin_stlist"])
@printf("     xcordat file:   %s\n", inpD["fin_xcordat"])
@printf("  travel time src:   %s\n", inpD["ttabsrc"])
@printf("     velocity mdl:   %s\n", inpD["fin_vzmdl"])
@printf("   map projection:   %s %s %.6f %.6f %.2f\n", 
    inpD["proj"], inpD["rellipse"], inpD["lon0"], inpD["lat0"], inpD["rotANG"])
@printf("\n")
@printf("****************** GROWCLUST Run Parameters *******************\n")
@printf("%56s %6.2f\n", " (min rxcor value for evpair similarity coeff.): rmin =", inpD["rmin"])
@printf("%56s %6.1f\n", " (max sta. dist for evpair similarity coeff.): delmax =", inpD["delmax"])
@printf("%56s %6.2f\n", " (max rms residual to join clusters): rmsmax =", inpD["rmsmax"])
@printf("%56s %6d\n" , " (num. bootstrap uncertainty iterations): nboot =", inpD["nboot"])
@printf("\n")
@printf("****************** Auxiliary Run Parameters *******************\n")
@printf("%56s %6.2f\n", " max catalog dist to join clusters: ", distmax)
@printf("%56s %6.2f\n", " max relocated dist to join clusters: ", distmax2)
@printf("%56s %6.2f\n", " max permitted horizontal cluster shifts: ", hshiftmax)
@printf("%56s %6.2f\n", " max permitted vertical cluster shifts: ", vshiftmax)
@printf("%56s %6.2f\n", " max median absolute residual to join clusters: ", rmedmax)  
@printf("***************************************************************\n")


# Run Summary with statistics
@printf( "\n")
@printf( "==================================================================\n")
@printf( "==================================================================\n")
@printf( "\n")
@printf( "********************  GROWCLUST Run Summary  *********************\n")
#@printf( "\n")
@printf( "%55s %10d\n", "Number of catalog events: ", nq)
@printf( "%55s %10d\n", "Number of relocated events: ", nreloc)
@printf( "%55s %10d\n", "Number of input event pairs: ", npair)
@printf( "%55s %10d\n", "Number of event pairs used: ", sum(qnpair)/2)
@printf( "%55s %10d\n", "Number of xcor data used (total, P+S): ", npp + nss)
@printf( "%55s %10d\n", "Number of xcor data used (P-phase): ", npp)
@printf( "%55s %10.4f\n", "RMS differential time residual (P-phase): ", rmsP)
@printf( "%55s %10.4f\n", "Mean (signed) differential time residual (P-phase): ", msresP)
@printf( "%55s %10d\n", "Number of xcor data used (S-phase): ", nss)
@printf( "%55s %10.4f\n", "RMS differential time residual (S-phase): ", rmsS)
@printf( "%55s %10.4f\n", "Mean (signed) differential time residual (S-phase): ", msresP)
@printf(  "\n")
@printf( "%55s %9d\n", "Number of clusters with >=   2 events: ", ntree2)
@printf( "%55s %9d\n", "Number of clusters with >=   5 events: ", ntree5)
@printf( "%55s %9d\n", "Number of clusters with >=  10 events: ", ntree10)
@printf( "%55s %9d\n", "Number of clusters with >=  20 events: ", ntree20)
@printf( "%55s %9d\n", "Number of clusters with >=  50 events: ", ntree50)
@printf( "%55s %9d\n", "Number of clusters with >= 100 events: ", ntree100)
   
### Done
println("\n\nDONE")
