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
const nupdate = 1000         # update progress every nupdate pairs
   
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
const ttabmode = "bystation" # travel time table mode ("bystation" or "byphase") - only for ray trace

# ------- Geodetic parameters -------------
const degkm = 111.1949266 # for simple degree / km conversion
const erad = 6371.0 # average earth radius in km, WGS-84
const mapproj = "tmerc" # projection for Proj4 ("aeqd", "lcc", "merc", "tmerc")
const rellipse = "WGS84" # reference ellipse for Proj4 (e.g. "WGS84")

### Read Input File

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
    inpD["itp_dz"] = inpD["tt_ddep"]
end


### Check Input parameters
input_ok = check_gcinp(inpD)
if input_ok
    println("Input parameters are ok!")
else
    println("ERROR: FIX INPUT PARAMETERS")
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
        boxwid, nit, irelonorm, vzmodel_type, mapproj)
if input_ok
    println("Auxiliary parameters are ok!")
else
    println("ERROR: FIX AUXILIARY PARAMETERS")
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
print("Travel-time table depths (min, max): ")
@printf("%.2f %.2f\n",inpD["tt_dep0"],inpD["tt_dep1"])
print("Travel-time table ranges (min, max): ")
@printf("%.2f %.2f\n",inpD["tt_del0"],inpD["tt_del1"])
println("Travel time grid directory: ",inpD["fdir_ttab"])
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
qlat0 = median(qdf.qlat)
qlon0 = median(qdf.qlon)
println("Median event location: $qlon0 $qlat0")

### Read Stations
print("\nReading station list")
@time sdf = read_stlist(inpD["fin_stlist"],inpD["stlist_fmt"])
min_selev = minimum(sdf.selev)
max_selev = maximum(sdf.selev)
mean_selev = mean(sdf.selev)
@printf("station elevation (min,mean,max): %.1fkm %.1fkm %.1fkm\n",
    min_selev, mean_selev, max_selev)
sta2elev = Dict(zip(sdf.sta,sdf.selev)) # map station to elevation

# ### Map Projection

# # setup projection
if mapproj in ["aeqd", "tmerc"]
    proj = Transformation("+proj=longlat +datum=$rellipse +no_defs",
        "+proj=$mapproj +datum=$rellipse +lat_0=$qlat0 +lon_0=$qlon0 +units=km")
elseif mapproj == "merc"
    proj = Transformation("+proj=longlat +datum=$rellipse +no_defs",
        "+proj=$mapproj +datum=$rellipse +lat_ts=$qlat0 +lon_0=$qlon0 +units=km")
elseif mapproj == "lcc"
    qlat1 = percentile(qdf.qlat,1.0)
    qlat2 = percentile(qdf.qlat,99.0)
    proj = Transformation("+proj=longlat +datum=$rellipse +no_defs",
    "+proj=$mapproj +datum=$rellipse +lat_0=$qlat0 +lon_0=$qlon0 +lat_1= $qlat1 +lat_2=$qlat2 +units=km")
else
    println("ERROR, map projection not defined! ", mapproj)
    exit()
end
iproj = inv(proj) # inverse transformation

# project station coordinates
sdf[!,:sX4] .= 0.0
sdf[!,:sY4] .= 0.0
for ii = 1:nrow(sdf)
    sdf[ii,[:sX4, :sY4]] = proj([sdf.slon[ii] sdf.slat[ii]])
end
println("\nProjected station list:")
show(sdf)
println()

# project event coordinates
qdf[!,:qX4] .= 0.0
qdf[!,:qY4] .= 0.0
for ii = 1:nrow(qdf)
    qdf[ii,[:qX4, :qY4]] = proj([qdf.qlon[ii] qdf.qlat[ii]])
end
println("\nProjected event list:")
show(qdf)
println()

### Read Xcor Data
println("\nReading xcor data") # in projected coordinates
@time xdf = read_xcordata_proj(inpD,qdf[!,[:qix,:qid,:qX4,:qY4]],sdf[!,[:sta,:sX4,:sY4]])
show(xdf)
println()

### Gather Unique Station/Phase combos
usta = unique(xdf[!,"sta"])
nstaU = length(usta)
println("$nstaU unique stations used:\n",usta,"\n")
xdf[!,:itab] .= Int16(0)
if (ttabmode == "byphase") & (inpD["ttabsrc"] == "trace") # one per phase
    xdf[!,:itab] .= convert.(Int16,xdf[!,:iphase]) # convert to Int16 for compatibility
    ntab = 2
else # one per station otherwise
    staIDX = Dict(zip(usta,1:length(usta)))
    for ii = 1:nrow(xdf) # table index: first nstaU are for P-waves, next are for S-waves
        xdf[ii,:itab] = staIDX[xdf[ii,:sta]] + nstaU*(xdf[ii,:iphase]-1)
    end
    ntab = 2*nstaU
end
println("Updated with table index:")
show(xdf[!,[:qix1,:qix2,:sta,:tdif,:rxcor,:iphase,:itab]])

###

###### Ray-Tracing Section ##########
if inpD["ttabsrc"] == "trace"

    ### Read in velocity model
    println("\nReading velocity model..")
    # read in
    z_s0, alpha_s0, beta_s0 = read_vzmodel(
        inpD["fin_vzmdl"],vpvs=inpD["vpvs_factor"])
    for ii in 1:length(z_s0) # print out
        @printf("%5.2fkm: %6.4f %6.4f\n",z_s0[ii],alpha_s0[ii],beta_s0[ii])
    end

    ### Find Moho depth in model, print results
    println("\nMoho depths:")
    imoho = find_moho(z_s0,alpha_s0)
    @printf("%.2f %.3f %.3f\n",
        z_s0[imoho], alpha_s0[imoho], beta_s0[imoho])
    println("Moho slownesses:")
    @printf("%.4f %.4f\n",1.0/alpha_s0[imoho], 1.0/beta_s0[imoho])

    # update minimum ray parameter
    if inpD["rayparam_min"] < 0.0
        inpD["plongcutP"] = 1.0/alpha_s0[imoho]
        inpD["plongcutS"] = 1.0/beta_s0[imoho]
    end

    ### Interpolate VZ Model to finer grid
    println("\nInterpolation:") # interpolate
    z_s, alpha_s, beta_s = interp_vzmodel(
        z_s0, alpha_s0, beta_s0; itp_dz=inpD["itp_dz"], ztop=-max_selev)
    for ii in 1:length(z_s) # print out
        @printf("%5.2fkm: %6.4f %6.4f\n",z_s[ii],alpha_s[ii],beta_s[ii])
    end

    ### Run Earth-flattening Codes
    z, alpha = eflatten(z_s, alpha_s, erad=erad)
    z, beta = eflatten(z_s, beta_s, erad=erad)

    ### Define slowness arrays
    npts = length(z_s)
    slow = zeros(Float64,npts,2)
    slow[:,1] .= 1.0./alpha
    slow[:,2] .= ifelse.(beta.>0.0,1.0./beta,1.0./alpha)

    ### Main Ray Tracing Loop
    println("\nRay tracing to assemble travel time tables:")

    # define grids
    qdeptab = collect(range(inpD["tt_dep0"],inpD["tt_dep1"],step=inpD["tt_ddep"]))
    sdeltab = collect(range(inpD["tt_del0"],inpD["tt_del1"],step=inpD["tt_ddel"]))

    # store all tables here
    ttLIST = [] # stations for P, then stations for S
    phases = ["P","S"] # make P and S tables
    plongcuts = [inpD["plongcutP"],inpD["plongcutS"]] # cutoffs

    # loop over phases
    println("STARTING RAY TRACING:")
    total_time = @elapsed for (iphase, phase) in enumerate(phases)

        # Print results
        println("Working on Phase $phase\n:")

        # loop over stations 
        for ista = 1:Int64(ntab/2)

            # get station depth (negative of elevation)
            if ttabmode =="bystation"
                sta = usta[ista] # station name
                zstart = -sta2elev[sta] # elevation
                println("Working on phase $phase station $sta: depth $zstart km")
                ttabfile = @sprintf("tt.%s.%sg",sta,lowercase(phase))
            else
                zstart = 0.0 # default is zero depth
                ttabfile = @sprintf("tt.%sg",lowercase(phase))
            end
        
            # Ray tracing: compute offset and travel time to different depths
            ptab, qdepxcor, qdeptcor, qdepucor, del2W, tt2W = trace_rays(
                iphase,z_s,z,slow,qdeptab,inpD["itp_dz"],zstart)
            
            # Compute slowness at station elevation
            isurf = findfirst(x->x>=zstart,z_s)
            if (isurf==1)|(z_s[isurf]==zstart)
                usurf=slow[isurf,iphase] # use slowness directly
            else                         # simple linear interpolation   
                usurf=slow[isurf-1,iphase] + (slow[isurf,iphase]-slow[isurf-1,iphase])*(
                    zstart-z_s[isurf-1])/(z_s[isurf]-z_s[isurf-1])
            end

            # Make table of first arrivals and take of angles
            TT, AA = first_arrivals(vzmodel_type, plongcuts[iphase], qdeptab, sdeltab, 
                                usurf, zstart, ptab, qdepxcor, qdeptcor, qdepucor, del2W, tt2W)

            # define interpolant object and add it to list
            iTT = make_smtrace_table(sdeltab,qdeptab,convert.(Float32,TT),shallowmode)
            push!(ttLIST,iTT)

            # optional output
            if inpD["fdir_ttab"] != "NONE"
                write_smtrace_table(inpD["fdir_ttab"] * ttabfile,inpD["fin_vzmdl"],iphase,vzmodel_type,
                            TT,qdeptab, sdeltab,ptab,zstart)
            end

        end # end loop over stations
        
    end
    
    @printf("\nElapsed seconds: %.2f",total_time)
    println()

    # finalize travel time tables
    const ttTABs = [deepcopy(xx) for xx in ttLIST] # static version
    const pTT = ttTABs[1] # p-wave example
    const sTT = ttTABs[Int64(ntab/2)+1] # s-wave example   

###### End of Ray Tracing Section ########

######## Pre-computed NonLinLoc Grids ########
elseif inpD["ttabsrc"] == "nllgrid"

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
            
            # update header
            grdparams["interp_mode"] = "linear"
            grdparams["xbounds"] = [inpD["tt_del0"],inpD["tt_del1"]]
            grdparams["zbounds"] = [inpD["tt_dep0"],inpD["tt_dep1"]]
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
    println("Travel time calculation mode not yet implemented: $ttsrc")
    println("Ending program.")
    exit()
end

##### Robustness checks for Tables/Grids #####

# find valid distances
println("\nChecking maximum allowable distance for ray tracing.")
maxdistTT = inpD["tt_del1"]
test_dists = collect(range(inpD["tt_del0"],inpD["tt_del1"],step=0.1))
dep1 = max(inpD["tt_dep0"],floor(minimum(qdf[!,"qdep"])))
dep2 = min(inpD["tt_dep1"],ceil(maximum(qdf[!,"qdep"])))
for test_depth in range(dep1,dep2,step=1.0)
    ttP = pTT(test_dists,test_depth)
    ttS = sTT(test_dists,test_depth)
    inan = (isnan.(ttP)) .| (isnan.(ttS))
    if sum(inan)>0
        global maxdistTT=min(maxdistTT,test_dists[inan][1]-1.0)
    end
end
println("Given travel time tables, maximum allowable distance: ",maxdistTT)

# define test distances, depth
println("Testing travel time interpolation (P and S waves):")
test_dists = collect(range(0.0,60.0,step=5.0))
npts = length(test_dists)
test_depth = 10.0
     
# calculate travel times
test_ttP = pTT(test_dists, test_depth)
test_ttS = sTT(test_dists, test_depth)

# print results
println(typeof(test_ttP[1]), " ", typeof(test_ttS[1]))
for ii in 1:npts     
@printf("Distance: %4.0f, Depth: %3.0f",
    test_dists[ii],test_depth)
@printf(" --> P and S travel times: %5.2f, %5.2f\n",
    test_ttP[ii],test_ttS[ii])
end

#### Validate Event Depths and Travel Time Tables

println("\nChecking event depths")

# print event depths
const min_qdep = minimum(qdf.qdep)
const max_qdep = maximum(qdf.qdep)
@printf("min and max event depth: %.3fkm %.3fkm\n",min_qdep,max_qdep)

# print table depths
@printf("min and max table depth: %.3fkm %.3fkm\n",inpD["tt_dep0"],inpD["tt_dep1"])

# implement warnings and checks
if (min_qdep < inpD["tt_dep0"])
    println("WARNING: min event depth < min table depth")
end
if (inpD["ttabsrc"]=="trace")
    if (min_qdep < z_s0[1]) # only checked for ray-trace
        println("WARNING: min event depth < min vzmodel depth")
        println("Results may be inaccurate near free-surface...")
        #exit() # allow this, but warn user (should be ok if depth is near 0)
    end
end
if (max_qdep > inpD["tt_dep1"]) # note tt_dep1 is >= vzmax
    println("ERROR: max event depth > max table / velocity model depth")
    exit()
end

println("Done.")


### Check xcor data
println("\nValidating xcor data...")
nbad = sum(abs.(xdf.tdif).>tdifmax)
if nbad > 0
    println("Error: bad input differential times, some larger than tdifmax=$tdifmax")
    println("Fix input file or adjust tdifmax parameter.")
    ibad = abs.(xdf.tdif).>tdifmax
    show(xdf[ibad,:])
    exit()
end
println("Max station distance: ",maximum(xdf.sdist))
nbad = sum(xdf.sdist.>inpD["tt_del1"])
if nbad > 0
    println("Error: bad input xcor data, stations further than travel time table allows")
    println("Fix input xcor file or adjust travel time table parameter at top of script.")
    ibad = xdf.sdist.>inpD["tt_del1"]
    show(xdf[ibad,:])
    exit()
elseif maxdistTT < inpD["tt_del1"]
    nbad = sum(xdf.sdist.>maxdistTT)
    if nbad > 0
        println("Error: bad input xcor data, stations further than allowed by ray tracing.")
        println("Given present configuration, rays cannot reach beyond a distance of $maxdistTT.")
        ibad = xdf.sdist.>maxdistTT
        show(xdf[ibad,:])
        exit()
    end 
end
println("Done.")

#### Done with robustness checks ###########


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
xdf00 = select(xdf,[:qix1,:qix2,:sX4,:sY4,:tdif,:itab,:igood])
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
    bpdf = combine(groupby(rxdf[!,Not([:sX4,:sY4,:itab,:tdif])],[:qix1,:qix2]),
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
    brXs, brYs, brdeps, brorgs, brcids, bnb = clustertree(
        bpdf.qix1, bpdf.qix2, bpdf.ix1, bpdf.ix2, 
        rxdf.tdif, rxdf.sX4, rxdf.sY4, rxdf.itab,
        qdf.qX4, qdf.qY4, qdf.qdep,
        ttTABs, nit, boxwid, irelonorm,
        inpD["rmsmax"],rmedmax,distmax,distmax2,
        hshiftmax,vshiftmax,torgdifmax,nupdate,maxlink)

    # inverse projection back to map coordinates
    brlons, brlats = zeros(nq), zeros(nq)
    for ii = 1:nq
        brlons[ii], brlats[ii] = iproj([brXs[ii] brYs[ii]])
    end
    
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
    @printf("Completed bootstrap iteration: %d/%d, wall clock = %.1fs.\n",ib,inpD["nboot"],wc)

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
resdf = select(xdf,[:qid1,:qid2,:tdif,:itab,:sX4,:sY4,:rxcor,:sdist])

# merge with event data, renaming columns to specify 1/2
resdf = innerjoin(resdf,rdf[!,[:enum,:evid,:rX,:rY,:rdep,:rtim]],on=:qid1=>:evid)
DataFrames.rename!(resdf,:enum=>:qnum1,:rX=>:qX1,:rY=>:qY1,:rdep=>:qZ1,:rtim=>:qtim1)
resdf = innerjoin(resdf,rdf[!,[:enum,:evid,:rX,:rY,:rdep,:rtim]],on=:qid2=>:evid)
DataFrames.rename!(resdf,:enum=>:qnum2,:rX=>:qX2,:rY=>:qY2,:rdep=>:qZ2,:rtim=>:qtim2)

# here, compute misfits for relocated event pairs and good differential times (for consistency with f90 codes)
resdf = resdf[(rcids[resdf.qnum1] .== rcids[resdf.qnum2]) .& ( # only event pairs in same cluster
    resdf[!,:sdist].<=inpD["delmax"]).&(resdf[!,:rxcor].>=inpD["rmin"]),:]

# compute individual source station distances
sdist1 = xydist(resdf[!,:qX1],resdf[!,:qY1],resdf[!,:sX4],resdf[!,:sY4])
sdist2 = xydist(resdf[!,:qX2],resdf[!,:qY2],resdf[!,:sX4],resdf[!,:sY4])

# compute predicted travel times
resdf[!,:pdif] .= 0.0
for ii = 1:nrow(resdf)
    resdf[ii,:pdif] = ttTABs[resdf[ii,:itab]](sdist2[ii],resdf[ii,:qZ2]) -
        ttTABs[resdf[ii,:itab]](sdist1[ii],resdf[ii,:qZ1]) +
        resdf[ii,:qtim2] - resdf[ii,:qtim1] # otime adjustment (add here or subtract from tdif)
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

### Compute bootstrap statistics ###

# pre-allocate: defaults are NaN for errors, 0 for nb arrays
boot_stdH, boot_madH = fill(NaN64,nq), fill(NaN64,nq)
boot_stdZ, boot_madZ = fill(NaN64,nq), fill(NaN64,nq)
boot_stdT, boot_madT = fill(NaN64,nq), fill(NaN64,nq)
boot_nbL, boot_nbM, boot_nbH = zeros(Int64,nq), zeros(Float64,nq), zeros(Int64,nq)

# loop over events
if inpD["nboot"] > 1 # need to have run bootstrapping


    println("\nComputing bootstrap statistics...")
    for ii = 1:nq

        # nbranch statistics
        boot_nbL[ii] = minimum(bnbM[ii,:])
        boot_nbM[ii] = mean(bnbM[ii,:])
        boot_nbH[ii] = maximum(bnbM[ii,:])
        
        # only if relocated
        if rdf[ii,:rnb]>1

            # standard errors
            xstd = degkm*cosd(rdf[ii,:rlat])*std(blonM[ii,:])
            ystd = degkm*std(blatM[ii,:])
            boot_stdH[ii]=sqrt(xstd^2+ystd^2)
            boot_stdZ[ii]=std(bdepM[ii,:])
            boot_stdT[ii]=std(borgM[ii,:])

            # MAD statistics (no division by 1/quantile(Normal(), 3/4) ≈ 1.4826)
            xmad = degkm*cosd(rdf[ii,:rlat])*mad(blonM[ii,:],normalize=false)
            ymad = degkm*mad(blatM[ii,:],normalize=false)
            boot_madH[ii]=sqrt(xmad^2+ymad^2)
            boot_madZ[ii]=mad(bdepM[ii,:],normalize=false)
            boot_madT[ii]=mad(borgM[ii,:],normalize=false)

        end
    end
end
        
#########################################################

### Write Output File: Catalog

println("\nWriting output catalog: ", inpD["fout_cat"])

# open output file
fcat = open(inpD["fout_cat"],"w")

# loop over events
for ii = 1:nq
    
    # print out origin time, relocated position, magnitude
    dateS = Dates.format(rdf[ii,:rot],"YYYY mm dd HH MM SS.sss")
    @printf(fcat,"%s %9d %9.5f %10.5f %7.3f %5.2f ",
        dateS,rdf[ii,:evid],rdf[ii,:rlat],rdf[ii,:rlon],
        rdf[ii,:rdep],qdf[ii,:qmag])
    
    # print out cluster number and fits
    @printf(fcat,"%7d %7d %7d %5d %5d %5d %5.2f %5.2f ",
        rdf[ii,:enum],rdf[ii,:rcid],rdf[ii,:rnb],qnpair[ii],
        qndiffP[ii],qndiffS[ii],qrmsP[ii],qrmsS[ii])
    
    # print out uncertanties and catalog locations
    @printf(fcat,"%7.3f %7.3f %7.3f %9.5f %10.5f %7.3f\n",
        boot_madH[ii],boot_madZ[ii],boot_madT[ii],
        qdf[ii,:qlat],qdf[ii,:qlon],qdf[ii,:qdep])
end

# close file
close(fcat)



### Write Output File: Cluster

if !(inpD["fout_clust"] in ["none","NONE"])

    println("\nWriting output clusters")

    # open output file
    fcc = open(inpD["fout_clust"],"w")

    # loop over all clusters to output (only makes sense for n>=2)
    for cc in findall(tnbranch .>= min(2,inpD["nbranch_min"]))
        
        # write cluster info
        @printf(fcc,"%8d %7d %9.5f %10.5f %7.3f %7.3f\n",
            cc,tnbranch[cc],tlats[cc],tlons[cc],tdeps[cc],torgs[cc])
        
        # write info for all events
        for ii in findall(rcids.==cc)
        
            # print out clustering, event info, mag, otime
            dateS = Dates.format(rdf[ii,:rot],"YYYY mm dd HH MM SS.sss")
            @printf(fcc,"%8d %8d %9d %5.2f %s ",cc,ii,rdf[ii,:evid],qdf[ii,:qmag],dateS)
            
            # print event location and position w/in cluster
            qdx = degkm*(rdf[ii,:rlon]-tlons[cc])*cosd(tlats[cc])
            qdy = degkm*(rdf[ii,:rlat]-tlats[cc])
            qdz = rdf[ii,:rdep]-tdeps[cc]
            @printf(fcc,"%9.5f %10.5f %7.3f %9.4f %9.4f %9.4f ",
                rdf[ii,:rlat],rdf[ii,:rlon],rdf[ii,:rdep],qdx,qdy,qdz)

            # print out uncertanties  and catalog locations
            @printf(fcc,"%7.3f %7.3f %7.3f %9.5f %10.5f %7.3f\n",
                boot_madH[ii],boot_madZ[ii],boot_madT[ii],
                qdf[ii,:qlat],qdf[ii,:qlon],qdf[ii,:qdep])
        end
        
    end


    # close file 
    close(fcc)

end

### Write Output File: Bootstrapping

# only if requested
if (inpD["nboot"] > 1)&(!(inpD["fout_boot"] in ["none","NONE"]))
    
    println("\nWriting output bootstrapping")
    
    # open file
    fbb = open(inpD["fout_boot"],"w")
    
    # file header
    @printf(fbb,"%8d %5d\n",nq,inpD["nboot"])
    
    
    # loop over events
    for ii = 1:nq
        
        # event header
        @printf(fbb,"%8d %9d %9.5f %10.5f %7.3f ",
            rdf[ii,:enum],rdf[ii,:evid],rdf[ii,:rlat],rdf[ii,:rlon],rdf[ii,:rdep])
        @printf(fbb,"%7d %7d %6.2f %7d %7d %9.5f %10.5f %7.3f ",
            rdf[ii,:rcid],rdf[ii,:rnb],boot_nbM[ii],boot_nbL[ii],boot_nbH[ii],
            qdf[ii,:qlat],qdf[ii,:qlon],qdf[ii,:qdep])
        @printf(fbb,"%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
            boot_madH[ii],boot_madZ[ii],boot_madT[ii],
            boot_stdH[ii],boot_stdZ[ii],boot_stdT[ii])
        
        # report each bootstrap (f90 version reports CID but this is dumb?)
        for ib = 1:inpD["nboot"]
            @printf(fbb,"%9.5f %10.5f %7.3f %7d\n",
                blatM[ii,ib],blonM[ii,ib],bdepM[ii,ib],bnbM[ii,ib])
        end
        
    end
    
    # close file
    close(fbb) 
    
end


### Write Output File: Log / Statistics

println("\nWriting run log")

# open log file
flog = open(inpD["fout_log"],"w")

# Run Parameters
@printf(flog,  "************************ Input files ************************\n")
@printf(flog,  "     control file:   %s\n", infile_ctl)
@printf(flog,  "       event list:   %s\n", inpD["fin_evlist"])
@printf(flog,  "     station list:   %s\n", inpD["fin_stlist"])
@printf(flog,  "     xcordat file:   %s\n", inpD["fin_xcordat"])
@printf(flog,  "     velo/grd mdl:   %s\n", inpD["fin_vzmdl"])
@printf(flog, "\n")
@printf(flog,  "************************ Output files *************************\n")
@printf(flog,  "     catalog file:   %s\n", inpD["fout_cat"])
@printf(flog,  "     cluster file:   %s\n", inpD["fout_clust"])
@printf(flog,  "         log file:   %s\n", inpD["fout_log"])
@printf(flog,  "   bootstrap file:   %s\n", inpD["fout_boot"])
@printf(flog, "\n")
@printf(flog,  "****************** GROWCLUST Run Parameters *******************\n")
@printf(flog, "%56s %6.2f\n", " (min rxcor value for evpair similarity coeff.): rmin =", inpD["rmin"])
@printf(flog, "%56s %6.1f\n", " (max sta. dist for evpair similarity coeff.): delmax =", inpD["delmax"])
@printf(flog, "%56s %6.2f\n", " (max rms residual to join clusters): rmsmax =", inpD["rmsmax"])
@printf(flog, "%56s %6d\n" , " (num. bootstrap uncertainty iterations): nboot =", inpD["nboot"])
@printf(flog, "\n")
@printf(flog,  "****************** Auxiliary Run Parameters *******************\n")
@printf(flog, "%56s %6.2f\n", " max catalog dist to join clusters: ", distmax)
@printf(flog, "%56s %6.2f\n", " max relocated dist to join clusters: ", distmax2)
@printf(flog, "%56s %6.2f\n", " max permitted horizontal cluster shifts: ", hshiftmax)
@printf(flog, "%56s %6.2f\n", " max permitted vertical cluster shifts: ", vshiftmax)
@printf(flog, "%56s %6.2f\n", " max median absolute residual to join clusters: ", rmedmax)  
@printf(flog,  "***************************************************************\n")


# Run Summary with statistics
@printf(flog, "\n")
@printf(flog, "==================================================================\n")
@printf(flog, "==================================================================\n")
@printf(flog, "\n")
@printf(flog, "********************  GROWCLUST Run Summary  *********************\n")
#@printf(flog, "\n")
@printf(flog, "%55s %10d\n", "Number of catalog events: ", nq)
@printf(flog, "%55s %10d\n", "Number of relocated events: ", nreloc)
@printf(flog, "%55s %10d\n", "Number of input event pairs: ", npair)
@printf(flog, "%55s %10d\n", "Number of event pairs used: ", sum(qnpair)/2)
@printf(flog, "%55s %10d\n", "Number of xcor data used (total, P+S): ", npp + nss)
@printf(flog, "%55s %10d\n", "Number of xcor data used (P-phase): ", npp)
@printf(flog, "%55s %10.4f\n", "RMS differential time residual (P-phase): ", rmsP)
@printf(flog, "%55s %10.4f\n", "Mean (signed) differential time residual (P-phase): ", msresP)
@printf(flog, "%55s %10d\n", "Number of xcor data used (S-phase): ", nss)
@printf(flog, "%55s %10.4f\n", "RMS differential time residual (S-phase): ", rmsS)
@printf(flog, "%55s %10.4f\n", "Mean (signed) differential time residual (S-phase): ", msresP)
@printf(flog,  "\n")
@printf(flog, "%55s %9d\n", "Number of clusters with >=   2 events: ", ntree2)
@printf(flog, "%55s %9d\n", "Number of clusters with >=   5 events: ", ntree5)
@printf(flog, "%55s %9d\n", "Number of clusters with >=  10 events: ", ntree10)
@printf(flog, "%55s %9d\n", "Number of clusters with >=  20 events: ", ntree20)
@printf(flog, "%55s %9d\n", "Number of clusters with >=  50 events: ", ntree50)
@printf(flog, "%55s %9d\n", "Number of clusters with >= 100 events: ", ntree100)
   
# close this
close(flog)
