### Script for Running GrowClust3D: Julia

###### Import packages ######

# External Packages
using Printf
using DataFrames
using Random
using Dates
using Proj: Transformation, inv
using Distributed
using SharedArrays

# GrowClust3D
using GrowClust3D

####### Define Algorithm Parameters (modify as necessary) #######

# ------- GrowClust algorithm parameters -------------------------------
const distmax = 8.0           # maximum catalog(input) distance to join clusters (km)
const distmax2 = 6.0          # maximum relocated distance to join clusters (km)
const hshiftmax = 2.0         # maximum permitted horizontal cluster shifts (km)
const vshiftmax = 2.0         # maximum permitted vertical cluster shifts (km)
const rmedmax = Float32(0.1) # maximum median absolute tdif residual to join clusters
const maxlink = 12            # use N best event pairs to relocate
const nbeststa = 24           # use N best xcorr values per event pair
const nupdate = 10000         # update progress every nupdate pairs
   
# ------- Relative Relocation subroutine parameters -------------
const boxwid = 3.                # initial "shrinking-box" width (km)
const nit = 15                   # number of iterations
const irelonorm = 1              # relocation norm (L1 norm=1, L2 norm=2, 3=robust L2)
const tdifmax = Float32(30.)     # maximum differential time value allowed (for error-checking on xcor data input)
const torgdifmax = Float32(10.0) # maximum origin time adjustment (for robustness)
   
# -------- Bootstrap resampling parameters -------------------
const iseed = 0 # random number seed
        
# ------- Velocity model parameters  -------------------------
const vzmodel_type = 1 # velocity model type: 1 = flat earth, (Z,Vp,Vs)
                 #                        or  2 = radial, (R,Vp,Vs): 
                 #               note: option 2 has not been extensively tested
const shallowmode = "flat" # option for how to treat shallow seismicity for 1D ray tracing
                       # flat treats negative depths as zero depth
                       # throw makes an error for negative depths
                       # linear does linear interpolation to negative depths
                       # reflect treats negative depths as equivalent to -depth

# ------- Geodetic parameters -------------
const degkm = 111.1949266 # for simple degree / km conversion
const erad = 6371.0       # average earth radius in km, WGS-84

###############################################################

############## Read and Organize Input Files ###################

# read input file
println("\nReading input file: ",ARGS[1])
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

### Print Finalized Input Parameters
print_input(inpD)

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
else
    plat1, plat2 = Nothing, Nothing # placeholders
end
if !(mlon01 < plon0 < mlon99)
    println("Caution! Projection origin outside bounds of seismicity.")
elseif !(mlat01 < plat0 < mlat99)
    println("Caution! Projection origin outside bounds of seismicity.")
end

### Read Stations
print("\nReading station list")
@time sdf = read_stlist(inpD["fin_stlist"],inpD["stlist_fmt"])
min_selev, max_selev, mean_selev = minimum(sdf.selev), maximum(sdf.selev), mean(sdf.selev)
@printf("station elevation (min,mean,max): %.3fkm %.3fkm %.3fkm\n",
    min_selev, mean_selev, max_selev)
sta2elev = Dict(zip(sdf.sta,sdf.selev)) # map station to elevation

##### Map Projections

## setup projection
proj = setup_projection(inpD,plon0,plat0,plat1,plat2)
iproj = inv(proj) # inverse transformation

# project station coordinates, including rotation
sX4, sY4 = lonlat2xypos(sdf.slon,sdf.slat,inpD["rotANG"],proj)
sdf[!,:sX4] .= sX4
sdf[!,:sY4] .= sY4

# project event coordinates
qX4, qY4 = lonlat2xypos(qdf.qlon,qdf.qlat,inpD["rotANG"],proj)
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

### Check for problems
nbad = sum(abs.(xdf.tdif).>tdifmax)
if nbad > 0 # validate differential times
    println("Error: bad input differential times, some larger than tdifmax=$tdifmax")
    println("Fix input file or adjust tdifmax parameter.")
    ibad = abs.(xdf.tdif).>tdifmax
    show(xdf[ibad,:])
    exit()
end

### Gather Unique Stations
usta = unique(xdf[!,"sta"])
nstaU = length(usta)
println("$nstaU unique stations used:\n",usta,"\n")

### Travel Time Table Index: One Per Station and Phase
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
@printf("Updated station list (only stations listed in xcorr file...)\n")
min_selev, max_selev, mean_selev = minimum(sdf.selev), maximum(sdf.selev), mean(sdf.selev)
@printf("min and max staZ: %.3fkm %.3fkm\n", min_selev, max_selev)
@printf("min and max staX: %.4f %.4f\n", minSX, maxSX)
@printf("min and max staY: %.4f %.4f\n", minSY, maxSY)
@printf("min and max staR: %.4f %.4f\n", minSR, maxSR)

####################################################

########### Compile Travel Time Tables #############

if inpD["ttabsrc"] == "trace" # 1D ray tracing
    ttLIST = make_trace1D_tables(inpD,maxSR,max_selev,usta,sta2elev,ntab,
        erad,vzmodel_type,shallowmode)
else  # 1D nonlinloc grids
    ttLIST = make_nllgrid_tables(inpD,maxSR,usta,shallowmode)
end

# finalize travel time tables
const ttTABs = [deepcopy(xx) for xx in ttLIST] # static version

### Validate Event Depths and Travel Time Tables

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
if (max_qdep > inpD["tt_zmax"]) # note tt_zmax is >= vzmax
    println("ERROR: max event depth > max table / velocity model depth")
    exit()
end
println("Done.")

#########################################################################

############# Main Clustering Loop: Including Bootstrapping ##############

# loading packages @everywhere
@everywhere using Printf
@everywhere using DataFrames
@everywhere using Random
@everywhere using Distributed
@everywhere using SharedArrays
@everywhere using GrowClust3D

# shared parameters - is there a better way?
@everywhere nboot = $(inpD["nboot"])
@everywhere nit = $nit
@everywhere boxwid = $boxwid
@everywhere degkm = $degkm
@everywhere irelonorm = $irelonorm
@everywhere rmsmax = $(inpD["rmsmax"])
@everywhere ngoodmin = $(inpD["ngoodmin"])
@everywhere rmedmax = $rmedmax
@everywhere distmax = $distmax
@everywhere distmax2 = $distmax2
@everywhere hshiftmax = $hshiftmax
@everywhere vshiftmax= $vshiftmax
@everywhere torgdifmax = $torgdifmax
@everywhere nupdate = $nupdate
@everywhere maxlink = $maxlink
@everywhere rmincut = $(inpD["rmincut"])
@everywhere nbeststa = $nbeststa 
@everywhere ttTABs = $ttTABs
@everywhere nq = $(nrow(qdf))
@everywhere iseed = $iseed
@everywhere tt_ndim = $(inpD["tt_ndim"])

# initial locations
@everywhere qX = $(qdf.qlat)
@everywhere qY = $(qdf.qlon)
@everywhere qZ = $(qdf.qdep)

# Setup bootstrapping matrices
bXM = SharedArray(repeat(qX,1,nboot+1))
bYM = SharedArray(repeat(qY,1,nboot+1))
bdepM = SharedArray(repeat(qZ,1,nboot+1))
borgM = SharedArray(zeros(Float32,(nq,nboot+1)))
bnbM = SharedArray(zeros(Int64,(nq,nboot+1)))
bcidM = SharedArray(repeat(Vector{Int32}(1:nq),1,nboot+1))
bnpairM = SharedArray(zeros(Int64,nboot+1))

# base xcor dataframe to sample from
xdf00 = select(xdf,[:qix1,:qix2,:sX4,:sY4,:tdif,:itab,:rxcor,:igood])
@everywhere xdf00 = $(xdf00)

# sampling vector
if nrow(xdf00)<typemax(Int32)
    @everywhere nxc = $(Int32(nrow(xdf00)))
    @everywhere ixc = Vector{Int32}(1:nxc)
else
    @everywhere nxc = $(nrow(xdf00))
    @everywhere ixc = Vector{Int64}(1:nxc)
end


#### loop over each bootstrapping iteration
println("\n\nStarting relocation estimates, workers=",workers())
@time @sync @distributed for ib in 0:nboot # need to call @sync to ensure all workers finish    
    
    # log thread id
    @printf("Starting bootstrap iteration: %d/%d\n",ib,nboot)
    
    # timer for this thread
    wc = @elapsed begin
    
    # bootstrapping: resample data before run
    Random.seed!(iseed + ib) # different for each run
    wc2 = @elapsed begin
    println("Initializing xcorr data and event pairs.")
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

    # calculate event pair similarity --> new version
    bpdf = combine(groupby(rxdf[!,Not([:itab,:tdif])],[:qix1,:qix2]),
         :rxcor => (x -> topNmeanpad(x,nbeststa,pad=rmincut/2.0)) => :rfactor,
         :ixx => first => :ix1,:ixx => last => :ix2,:igood => sum => :ngood)
    
    # sort pairs (note, resampled pairs may not have ngoodmin tdifs)
    if ib > 0
        bpdf = bpdf[bpdf.ngood.>=ngoodmin-2,[:qix1,:qix2,:rfactor,:ix1,:ix2]] # for robustness
    else
        select!(bpdf,[:qix1,:qix2,:rfactor,:ix1,:ix2])
    end
    sort!(bpdf,:rfactor,rev=true) # so best pairs first
    end # ends elapsed time for setup
    println("Done with initialization, elapsed time = $wc2")
    
    # run clustering
    if tt_ndim < 3 # 2D travel time table
        brXs, brYs, brdeps, brorgs, brcids, bnb = clustertree(
            bpdf.qix1, bpdf.qix2, bpdf.ix1, bpdf.ix2, rxdf.tdif, rxdf.sX4, rxdf.sY4, 
            rxdf.itab, qdf.qX4, qdf.qY4, qdf.qdep,ttTABs, nit, boxwid, irelonorm,
            rmsmax,rmedmax,distmax,distmax2,hshiftmax,vshiftmax,torgdifmax,nupdate,maxlink)
    else          # 3D travel time table
        brXs, brYs, brdeps, brorgs, brcids, bnb = clustertree_3D(
            bpdf.qix1, bpdf.qix2, bpdf.ix1, bpdf.ix2, rxdf.tdif, rxdf.itab,
            qdf.qX4, qdf.qY4, qdf.qdep, ttTABs, nit, boxwid, irelonorm, rmsmax,
            rmedmax,distmax,distmax2,hshiftmax,vshiftmax,torgdifmax,nupdate,maxlink)
    end

    # save output to Shared Array
    if ib > 0
        bXM[:,ib] .= brXs
        bYM[:,ib] .= brYs
        bdepM[:,ib] .= brdeps
        borgM[:,ib] .= brorgs
        bnbM[:,ib] .= bnb
        bcidM[:,ib] .= brcids
        bnpairM[ib] = nrow(bpdf)
    else
        bXM[:,nboot+1] .= brXs
        bYM[:,nboot+1] .= brYs
        bdepM[:,nboot+1] .= brdeps
        borgM[:,nboot+1] .= brorgs
        bnbM[:,nboot+1] .= bnb
        bcidM[:,nboot+1] .= brcids
        bnpairM[nboot+1] = nrow(bpdf)
    end
        
    # completion
    end # ends the wall clock
    @printf("Completed bootstrap iteration: %d/%d, wall clock = %.1fs.\n",ib,nboot,wc)

end

### Invert Map Projection
blatM, blonM = zeros(Float64,nq,nboot+1), zeros(Float64,nq,nboot+1)
for jj=1:nboot+1
    blonM[:,jj], blatM[:,jj] = xypos2latlon(bXM[:,jj],bYM[:,jj],inpD["rotANG"],iproj)
end

### Extract relocated event-based output arrays
revids = qdf[:,:qid]
rXs, rYs = bXM[:,nboot+1], bYM[:,nboot+1]
rlats, rlons, rdeps = blatM[:,nboot+1], blonM[:,nboot+1], bdepM[:,nboot+1]
rorgs, rcids, npair = borgM[:,nboot+1], bcidM[:,nboot+1], bnpairM[nboot+1]

################################################################

### Finalize Clustering Trees and Relocated Catalog

# compile clustering tree
rnbranch,rcids,tnbranch,tlats,tlons,tdeps,torgs = make_clustertree(
    rlons,rlats,rdeps,rorgs,rcids,qdf.qix)

# initialize relocated dataframe
println("\nFinalizing relocated dataset:")
rdf = DataFrame("enum"=>qdf.qix,"evid"=>revids,"mag"=>qdf.qmag,
    "rlat"=>rlats,"rlon"=>rlons,"rdep"=>rdeps,"rtim"=>rorgs,
    "rcid"=>rcids,"rnb"=>rnbranch,"rX"=>rXs,"rY"=>rYs,
    "qlat"=>qdf.qlat,"qlon"=>qdf.qlon,"qdep"=>qdf.qdep)

# compute relocated
nreloc = sum(rdf.rnb .> 1)
println("Relocated: ",nreloc)

# adjusted otime
rdf[!,:rot] = qdf[!,:qotime] .+ Nanosecond.(round.(rorgs*1.0e9))
show(rdf)
println()

### Compute Misfits - w/ otime adjustment ###
npp, nss, rmsP, rmsS, msresP, msresS, qrmsP, qrmsS,
     qndiffP, qndiffS, qnpair = compute_misfits(inpD,xdf,rdf,ttTABs)

### Compute bootstrap statistics ###
boot_madH, boot_madZ, boot_madT, boot_stdH, boot_stdZ, boot_stdT,
    boot_nbL, boot_nbM, boot_nbH = compute_bootstats(inpD,degkm,rdf,bnbM,blonM,blatM,bdepM,borgM)
        
#########################################################

### Write Output File: Catalog (if requested)
if !(inpD["fout_cat"] in ["none","None", "NONE"])
    write_cat(inpD,rdf,qnpair,qndiffP,qndiffS,
        qrmsP,qrmsS,boot_madH,boot_madZ,boot_madT)
end

### Write Output File: Cluster (if requested)
if !(inpD["fout_clust"] in ["none","None", "NONE"])
    write_clust(inpD,rdf,tnbranch,tlats,tlons,tdeps,torgs,
        boot_madH, boot_madT, boot_madZ)
end

### Write Output File: Bootstrapping (if requested)
if (inpD["nboot"] > 1)&(!(inpD["fout_boot"] in ["none","None", "NONE"]))
    write_boot(inpD,rdf,boot_madH,boot_madT,boot_madZ,
        boot_nbH,boot_nbL,boot_nbM,boot_stdH,boot_stdT,boot_stdZ,
        blatM,blonM,bdepM,bnbM)
end

### Write Output File: Log / Statistics (or print to screen)
write_log(inpD,[infile_ctl, distmax, distmax2, hshiftmax, vshiftmax, rmedmax],
    [nq, nreloc, npair, qnpair, npp, nss, rmsP, rmsS, msresP, msresS], tnbranch)

### Report completion
@printf("\nCompleted task: run_growclust3D-MP.jl\n")