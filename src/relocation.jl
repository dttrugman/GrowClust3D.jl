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

### Robust Mean Based on Huber norm

function robomean(xx::Vector{Float64},xgap::Float64,nit::Int64)
    
    # get npts
    ndat = length(xx)
    
    # set initial weight
    xw = ones(Float64,ndat)
    
    # iterations loop
    xmean = 0.0
    for it = 1:nit
        
        # compute weighted mean
        xmean = sum(xw.*xx)/sum(xw)
        
        # compute weights for next iteration
        if it < nit # distance is xgap + |xx-xmean|, w = 1/d
            xw .= 1.0 ./ (xgap.+abs.(xx.-xmean))
        end
        
    end
    
    # compute final weighted misfit
    fit2 = sum(xw .* (xx.-xmean).^2)
    
    # return
    return xmean, fit2 
end

function robomean(xx::Vector{Float32},xgap::Float32,nit::Int64)
    
    # get npts
    ndat = length(xx)
    
    # set initial weight
    xw = ones(Float32,ndat)
    
    # iterations loop
    xmean = Float32(0.0)
    for it = 1:nit
        
        # compute weighted mean
        xmean = sum(xw.*xx)/sum(xw)
        
        # compute weights for next iteration
        if it < nit # distance is xgap + |xx-xmean|, w = 1/d
            xw .= Float32(1.0) ./ (xgap.+abs.(xx.-xmean))
        end
        
    end
    
    # compute final weighted misfit
    fit2 = sum(xw .* (xx.-xmean).^2)
    
    # return
    return xmean, fit2 
end

### Define rms function
function evalrms(data)
    return sqrt.(mean(data.^2))
end


# DIFCLUST performs relative relocation of two clusters of events
# (relative to the centroid of the cluster pair) using the the npr "best" 
# event pairs linking the clusters. (npr is not input explicitly,
# as the npick differential travel times, etc., include all observations 
# across these linking pairs). 
# This is analogous to DIFLOC, which does the same thing for two events, 
# relative to their centroid... As in DIFLOC, the method uses an iterative
# ("shrinking-box") grid search approach to obtain the best (L1) relative locations.
#----------
# Three versions: difclust1, difclust2, difclust3 --> for L1, L2, and Robomean norm.

# Inputs: qlat0  =  reference center point latitude
#         qlon0  =  reference center point longitude
#         qdep0  =  reference center point depth (km)
#         tdif   =  array (len=npick) of dif times, t2-t1 (s)
#         iph    =  array (len=npick) with phase index numbers (1 to 10) for tt data
#         slat   =  array (len=npick) with station latitudes
#         slon   =  array (len=npick) with station longitudes
#         qlat1  =  array (len=npick) of events in cluster1 latitude offsets from centroid
#         qlon1  =  array (len=npick) of events in cluster1 longitude offsets
#         qdep1  =  array (len=npick) of events in cluster1 depth offsets
#         qorg1  =  array (len=npick) of events in cluster1 time offsets
#         qlat2  =  array (len=npick) of events in cluster2 latitude offsets from centroid
#         qlon2  =  array (len=npick) of events in cluster2 longitude offsets
#         qdep2  =  array (len=npick) of events in cluster2 depth offsets
#         qorg2  =  array (len=npick) of events in cluster2 time offsets
#         ttTABs =  travel time tables for P and S phases
#         boxwid =  starting box width (km)
#         nit    =  number of iterations to perform
#         degkm  =  degrees to km conversion factor
# Returns: clat1  =  best-fitting latitude for first cluster centroid
#          clon1  =  best-fitting longitude for first cluster centroid
#          cdep1  =  best-fitting depth (km) of first cluster centroid
#          clat2  =  best-fitting latitude for second cluster centroid
#          clon2  =  best-fitting longitude for second cluster centroid
#          cdep2  =  best-fitting depth (km) of second cluster centroid
#          cdist  =  cluster separation distance (km)
#          torgdif=  origin time difference, i.e., t2-t1 median residual (>0 when 2 is later than 1)
#          resid  =  array (len=npick) of residuals (s) between observed tdif (tt) and predicted
#          rms    =  rms residual ( sqrt( sum(resid**2)/npick) )
#          rmed   =  median absolute value residual ( median( abs(resid) ) )
#          resol  =  nominal resolution (m) of final box

##### difclust1: L1 residual norm
function difclust1(qlat0::Float64,qlon0::Float64,qdep0::Float64,
    tdif::Vector{Float64},iph::Vector{Int8},
    slat::Vector{Float64},slon::Vector{Float64},
    qlat1::Vector{Float64},qlon1::Vector{Float64},
    qdep1::Vector{Float64},qorg1::Vector{Float64},
    qlat2::Vector{Float64},qlon2::Vector{Float64},
    qdep2::Vector{Float64},qorg2::Vector{Float64},
    ttTABs,boxwid::Float64,nit::Int64,degkm::Float64)

# initialize some variables
flatbest1, flonbest1, fdepbest1 = 0.0, 0.0, 0.0
flatbest2, flonbest2, fdepbest2 = 0.0, 0.0, 0.0
torgdif = 0.0

# extract npick
npick = length(tdif)
resid = zeros(npick)

# initialize box
dlat0 = 0.0
dlon0 = 0.0
ddep0 = 0.0
dlat = 0.5*boxwid/degkm
cosqlat = cosd(qlat0)
dlon = dlat/cosqlat

# careful near surface (added 07/2018)
zboxwid = boxwid
if (qdep0 < 0.0)
    zboxwid = min(1.0, boxwid)
elseif (qdep0 < boxwid/2.0)
    zboxwid = max(1.0,2.0*qdep0)
end
ddep = 0.5*zboxwid

# start iterations
for it=1:nit
    
    # grid search over box for best fit
    fitbest = 1.0e20
    tbest = 0.0
    
    # get trial locations: f1, f2 (3x3 grid for this box iteration)
    for iy = -1.0:1.0
        flat1 = qlat0 + dlat0 + dlat*iy
        flat2 = qlat0 - dlat0 - dlat*iy
        for ix = -1.0:1.0
            flon1 = qlon0 + dlon0 + dlon*ix
            flon2 = qlon0 - dlon0 - dlon*ix
            for iz = -1.0:1.0
                fdep1 = qdep0 + ddep0 + ddep*iz
                fdep2 = qdep0 - ddep0 - ddep*iz
                
                # compute predicted travel time and residuals w/observed
                sdist1 = map_distance(qlat1.+flat1,qlon1.+flon1,slat,slon)
                sdist2 = map_distance(qlat2.+flat2,qlon2.+flon2,slat,slon)                
                @inbounds for ii=1:npick
                    tt1 = ttTABs[iph[ii]](sdist1[ii],qdep1[ii]+fdep1)
                    tt2 = ttTABs[iph[ii]](sdist2[ii],qdep2[ii]+fdep2)
                    pdif = (tt2 + qorg2[ii]) - (tt1 + qorg1[ii])
                    resid[ii] = tdif[ii] - pdif # accounts for otime adjustment
                end
                
                # compute fit: L1 norm
                residval = median(resid)
                fit = sum(abs.(resid.-residval))
                
                # update best fit
                if fit < fitbest
                    fitbest = fit
                    flatbest1 = flat1
                    flonbest1 = flon1
                    fdepbest1 = fdep1
                    flatbest2 = flat2
                    flonbest2 = flon2
                    fdepbest2 = fdep2
                    tbest = residval
                end
                
            end
        end
    end # end of grid search
    
    # update best position
    dlat0 = flatbest1 - qlat0
    dlon0 = flonbest1 - qlon0
    ddep0 = fdepbest1 - qdep0
    torgdif = tbest # otime difference: residual between observed and predicted tdif

    # shrink box by 2/3 each iteration
    dlat *= 2.0/3.0      
    dlon *= 2.0/3.0
    ddep *= 2.0/3.0
    
end # end loop over iterations

# output final locations
clat1 = flatbest1
clon1 = flonbest1
cdep1 = fdepbest1
clat2 = flatbest2
clon2 = flonbest2
cdep2 = fdepbest2
resol = (dlat/(2.0/3.0))*degkm

# compute distance between cluster centroids
cdist = sqrt((cdep2-cdep1)^2 + (map_distance(clat1,clon1,clat2,clon2))^2)
                
# compute residual between observed and predicted travel time
sdist1 = map_distance(clat1.+qlat1,clon1.+qlon1,slat,slon)
sdist2 = map_distance(clat2.+qlat2,clon2.+qlon2,slat,slon)
@inbounds for ii=1:npick
    tt1 = ttTABs[iph[ii]](sdist1[ii],cdep1+qdep1[ii])
    tt2 = ttTABs[iph[ii]](sdist2[ii],cdep2+qdep2[ii])
    pdif = (tt2 + qorg2[ii]) - (tt1 + qorg1[ii])
    resid[ii] = tdif[ii] - pdif - torgdif
end

# rms residual and median absolute residual
rms = sqrt.(mean(resid.^2))
rmed = mad(resid,normalize=false) # default is to equate to sdev

# return results
return clat1, clon1, cdep1, clat2, clon2, cdep2, cdist, torgdif, resid, rms, rmed, resol

end

##### difclust2: L2 residual norm
function difclust2(qlat0::Float64,qlon0::Float64,qdep0::Float64,
    tdif::Vector{Float64},iph::Vector{Int8},
    slat::Vector{Float64},slon::Vector{Float64},
    qlat1::Vector{Float64},qlon1::Vector{Float64},
    qdep1::Vector{Float64},qorg1::Vector{Float64},
    qlat2::Vector{Float64},qlon2::Vector{Float64},
    qdep2::Vector{Float64},qorg2::Vector{Float64},
    ttTABs,boxwid::Float64,nit::Int64,degkm::Float64)

# initialize some variables
flatbest1, flonbest1, fdepbest1 = 0.0, 0.0, 0.0
flatbest2, flonbest2, fdepbest2 = 0.0, 0.0, 0.0
torgdif = 0.0

# extract npick
npick = length(tdif)
resid = zeros(npick)

# initialize box
dlat0 = 0.0
dlon0 = 0.0
ddep0 = 0.0
dlat = 0.5*boxwid/degkm
cosqlat = cosd(qlat0)
dlon = dlat/cosqlat

# careful near surface (added 07/2018)
zboxwid = boxwid
if (qdep0 < 0.0)
    zboxwid = min(1.0, boxwid)
elseif (qdep0 < boxwid/2.0)
    zboxwid = max(1.0,2.0*qdep0)
end
ddep = 0.5*zboxwid

# start iterations
for it=1:nit
    
    # grid search over box for best fit
    fitbest = 1.0e20
    tbest = 0.0
    
    # get trial locations: f1, f2 (3x3 grid for this box iteration)
    for iy = -1.0:1.0
        flat1 = qlat0 + dlat0 + dlat*iy
        flat2 = qlat0 - dlat0 - dlat*iy
        for ix = -1.0:1.0
            flon1 = qlon0 + dlon0 + dlon*ix
            flon2 = qlon0 - dlon0 - dlon*ix
            for iz = -1.0:1.0
                fdep1 = qdep0 + ddep0 + ddep*iz
                fdep2 = qdep0 - ddep0 - ddep*iz
                
                # compute predicted travel time and residuals w/observed
                sdist1 = map_distance(qlat1.+flat1,qlon1.+flon1,slat,slon)
                sdist2 = map_distance(qlat2.+flat2,qlon2.+flon2,slat,slon)                
                @inbounds for ii=1:npick
                    tt1 = ttTABs[iph[ii]](sdist1[ii],qdep1[ii]+fdep1)
                    tt2 = ttTABs[iph[ii]](sdist2[ii],qdep2[ii]+fdep2)
                    pdif = (tt2 + qorg2[ii]) - (tt1 + qorg1[ii])
                    resid[ii] = tdif[ii] - pdif # accounts for otime adjustment
                end
                
                # compute fit with L2 norm
                residval = mean(resid)
                fit = sum((resid.-residval).^2)
                
                # update best fit
                if fit < fitbest
                    fitbest = fit
                    flatbest1 = flat1
                    flonbest1 = flon1
                    fdepbest1 = fdep1
                    flatbest2 = flat2
                    flonbest2 = flon2
                    fdepbest2 = fdep2
                    tbest = residval
                end
                
            end
        end
    end # end of grid search
    
    # update best position
    dlat0 = flatbest1 - qlat0
    dlon0 = flonbest1 - qlon0
    ddep0 = fdepbest1 - qdep0
    torgdif = tbest # otime difference: residual between observed and predicted tdif

    # shrink box by 2/3 each iteration
    dlat *= 2.0/3.0      
    dlon *= 2.0/3.0
    ddep *= 2.0/3.0
    
end # end loop over iterations

# output final locations
clat1 = flatbest1
clon1 = flonbest1
cdep1 = fdepbest1
clat2 = flatbest2
clon2 = flonbest2
cdep2 = fdepbest2
resol = (dlat/(2.0/3.0))*degkm

# compute distance between cluster centroids
cdist = sqrt((cdep2-cdep1)^2 + (map_distance(clat1,clon1,clat2,clon2))^2)
                
# compute residual between observed and predicted travel time
sdist1 = map_distance(clat1.+qlat1,clon1.+qlon1,slat,slon)
sdist2 = map_distance(clat2.+qlat2,clon2.+qlon2,slat,slon)
@inbounds for ii=1:npick
    tt1 = ttTABs[iph[ii]](sdist1[ii],cdep1+qdep1[ii])
    tt2 = ttTABs[iph[ii]](sdist2[ii],cdep2+qdep2[ii])
    pdif = (tt2 + qorg2[ii]) - (tt1 + qorg1[ii])
    resid[ii] = tdif[ii] - pdif - torgdif
end

# rms residual and median absolute residual
rms = sqrt.(mean(resid.^2))
rmed = mad(resid,normalize=false) # default is to equate to sdev

# return results
return clat1, clon1, cdep1, clat2, clon2, cdep2, cdist, torgdif, resid, rms, rmed, resol

end

##### difclust3: L3 (robust mean) residual norm
function difclust3(qlat0::Float64,qlon0::Float64,qdep0::Float64,
    tdif::Vector{Float64},iph::Vector{Int8},
    slat::Vector{Float64},slon::Vector{Float64},
    qlat1::Vector{Float64},qlon1::Vector{Float64},
    qdep1::Vector{Float64},qorg1::Vector{Float64},
    qlat2::Vector{Float64},qlon2::Vector{Float64},
    qdep2::Vector{Float64},qorg2::Vector{Float64},
    ttTABs,boxwid::Float64,nit::Int64,degkm::Float64)

# initialize some variables
flatbest1, flonbest1, fdepbest1 = 0.0, 0.0, 0.0
flatbest2, flonbest2, fdepbest2 = 0.0, 0.0, 0.0
torgdif = 0.0

# extract npick
npick = length(tdif)
resid = zeros(npick)

# initialize box
dlat0 = 0.0
dlon0 = 0.0
ddep0 = 0.0
dlat = 0.5*boxwid/degkm
cosqlat = cosd(qlat0)
dlon = dlat/cosqlat

# careful near surface (added 07/2018)
zboxwid = boxwid
if (qdep0 < 0.0)
    zboxwid = min(1.0, boxwid)
elseif (qdep0 < boxwid/2.0)
    zboxwid = max(1.0,2.0*qdep0)
end
ddep = 0.5*zboxwid

# start iterations
for it=1:nit
    
    # grid search over box for best fit
    fitbest = 1.0e20
    tbest = 0.0
    
    # get trial locations: f1, f2 (3x3 grid for this box iteration)
    for iy = -1.0:1.0
        flat1 = qlat0 + dlat0 + dlat*iy
        flat2 = qlat0 - dlat0 - dlat*iy
        for ix = -1.0:1.0
            flon1 = qlon0 + dlon0 + dlon*ix
            flon2 = qlon0 - dlon0 - dlon*ix
            for iz = -1.0:1.0
                fdep1 = qdep0 + ddep0 + ddep*iz
                fdep2 = qdep0 - ddep0 - ddep*iz
                
                # compute predicted travel time and residuals w/observed
                sdist1 = map_distance(qlat1.+flat1,qlon1.+flon1,slat,slon)
                sdist2 = map_distance(qlat2.+flat2,qlon2.+flon2,slat,slon)                
                @inbounds for ii=1:npick
                    tt1 = ttTABs[iph[ii]](sdist1[ii],qdep1[ii]+fdep1)
                    tt2 = ttTABs[iph[ii]](sdist2[ii],qdep2[ii]+fdep2)
                    pdif = (tt2 + qorg2[ii]) - (tt1 + qorg1[ii])
                    resid[ii] = tdif[ii] - pdif # accounts for otime adjustment
                end
                
                # compute fit (depends on norm)
                residval, fit = robomean(resid,0.1,10)
                
                # update best fit
                if fit < fitbest
                    fitbest = fit
                    flatbest1 = flat1
                    flonbest1 = flon1
                    fdepbest1 = fdep1
                    flatbest2 = flat2
                    flonbest2 = flon2
                    fdepbest2 = fdep2
                    tbest = residval
                end
                
            end
        end
    end # end of grid search
    
    # update best position
    dlat0 = flatbest1 - qlat0
    dlon0 = flonbest1 - qlon0
    ddep0 = fdepbest1 - qdep0
    torgdif = tbest # otime difference: residual between observed and predicted tdif

    # shrink box by 2/3 each iteration
    dlat *= 2.0/3.0      
    dlon *= 2.0/3.0
    ddep *= 2.0/3.0
    
end # end loop over iterations

# output final locations
clat1 = flatbest1
clon1 = flonbest1
cdep1 = fdepbest1
clat2 = flatbest2
clon2 = flonbest2
cdep2 = fdepbest2
resol = (dlat/(2.0/3.0))*degkm

# compute distance between cluster centroids
cdist = sqrt((cdep2-cdep1)^2 + (map_distance(clat1,clon1,clat2,clon2))^2)
                
# compute residual between observed and predicted travel time
sdist1 = map_distance(clat1.+qlat1,clon1.+qlon1,slat,slon)
sdist2 = map_distance(clat2.+qlat2,clon2.+qlon2,slat,slon)
@inbounds for ii=1:npick
    tt1 = ttTABs[iph[ii]](sdist1[ii],cdep1+qdep1[ii])
    tt2 = ttTABs[iph[ii]](sdist2[ii],cdep2+qdep2[ii])
    pdif = (tt2 + qorg2[ii]) - (tt1 + qorg1[ii])
    resid[ii] = tdif[ii] - pdif - torgdif
end

# rms residual and median absolute residual
rms = sqrt.(mean(resid.^2))
rmed = mad(resid,normalize=false) # default is to equate to sdev

# return results
return clat1, clon1, cdep1, clat2, clon2, cdep2, cdist, torgdif, resid, rms, rmed, resol

end

#####################################

### CLUSTERTREE: Runs the the main growclust algorithm, relocated events using xcorr data / clustering.
#
# Inputs:
# - pqix1, pqix2: vectors of serial event numbers for each event pair
# - ixx1, ixx2: start and stop indices in xcorr arrays for each pair
# - tdif, slat, slon, iphase: differential travel time, station lat, station lon, phase for each xcorr data
# - qlats, qlons, qdeps: initial event locations
# - ttTABS: P and S travel time interpolation objects
# - nit, boxwid, degkm, irelonorm: difclust parameters
# - rmsmax, rmedmax, distmax, distmax2, hshiftmax, vshifmax, torgdifmax: cluster merge parameters
# - nupdate: counter for iteration update
# - maxlink: maximum number of event pairs used to link clusters
#
# Returns:
# - brlats, brlons, brdeps: relocated event locations
# - brcids, bnb: cluster ids and number of events in each cluster 
#
# function with i32 ixx arrays
function clustertree(pqix1::Vector{Int32},pqix2::Vector{Int32},ixx1::Vector{Int32},ixx2::Vector{Int32},
    tdif::Vector{Float64},slat::Vector{Float64},slon::Vector{Float64},iphase::Vector{Int8},
    qlats::Vector{Float64},qlons::Vector{Float64},qdeps::Vector{Float64},ttTABs,
    nit::Int64,boxwid::Float64,degkm::Float64,irelonorm::Int64,rmsmax::Float64,rmedmax::Float64,
    distmax::Float64,distmax2::Float64,hshiftmax::Float64,vshiftmax::Float64,torgdifmax::Float64,
    nupdate::Int64,maxlink::Int64)

    # setup parameters
    cdepmin = min(0.0,minimum(qdeps))
    hshiftmaxD2 = (hshiftmax/degkm)^2 # in squared degrees
    distmax22 = distmax^2 # in squared km

    # array sizes
    nq = length(qlats)
    bnpair = length(pqix1)

    # initialize relocated event arrays
    brlats = copy(qlats)
    brlons = copy(qlons)
    brdeps = copy(qdeps)
    brorgs = zeros(nq)
    brcids = Vector{Int32}(1:nq)

    # initialize clustering tree arrays
    btlats = copy(brlats)
    btlons = copy(brlons)
    btdeps = copy(brdeps)
    btorgs = zeros(nq)
    btnbranch = ones(Int32,nq)

    # dictionary with cid => event indices (note cid=qix on initialization)
    cid2qixD = Dict(qix=>[qix] for qix in brcids)

    # dictionary with cid => event pair indices
    cid2pairD = Dict( qix => Vector{Int32}(findall(
    (pqix1.==qix).|(pqix2 .== qix))) for qix in brcids) # cid=qix on initialization


    # loop over event pairs
    @inbounds for ip in Vector{Int32}(1:bnpair)

    # Progress
    if (mod(ip,nupdate)==0)
        println("Thread: ", Threads.threadid(), "--> Working on sorted pair: $ip/$bnpair")
    end

    # get event pair information
    qix1, qix2 = pqix1[ip], pqix2[ip]
    qc1, qc2 = brcids[qix1], brcids[qix2]
    if qc1 == qc2; continue; end # skip if in same cluster
    nb1, nb2 = btnbranch[qc1], btnbranch[qc2]

    # check to see if clusters are too far apart
    cdist22 = (btdeps[qc2]-btdeps[qc1])^2 + map_distance(
    btlats[qc1],btlons[qc1],btlats[qc2],btlons[qc2])^2
    if cdist22 > distmax22
        continue
    end

    # find event ids in each cluster
    c1qixs = cid2qixD[qc1] 
    c2qixs = cid2qixD[qc2]

    # only one link
    if (nb1==1)&(nb2==1)

        # only one link
        nlink = 1
        linx = [ip]

    # find all links between clusters
    else

        # # event pairs linking clusters
        # #   - pairs above are either processed and in same cluster
        # #     or caused issues when trying to link
        cpix1, cpix2 = cid2pairD[qc1],cid2pairD[qc2]
        @views linx = intersect(cpix1[cpix1.>=ip],cpix2[cpix2.>=ip])

        # keeping only best linx
        nlink = length(linx)
        if nlink > maxlink
            linx = partialsort(linx,1:maxlink)
        end 

    end

    # count number of picks
    npick = sum([ixx2[jj] - ixx1[jj] + 1 for jj in linx])

    # extract locations relative to centroid
    dqlat1, dqlon1 = zeros(npick),zeros(npick)
    dqdep1, dqorg1 = zeros(npick),zeros(npick)
    dqlat2, dqlon2 = zeros(npick),zeros(npick)
    dqdep2, dqorg2 = zeros(npick),zeros(npick)
    phase21, slat21, slon21 = zeros(Int8,npick), zeros(npick), zeros(npick)
    tdif21 = zeros(npick)
    ix1 = 1 # start index for event pair in the arrays above
    @inbounds for jj in linx
        pix1, pix2, jx1, jx2 = pqix1[jj],pqix2[jj],ixx1[jj],ixx2[jj]
        ix2 = ix1 + (jx2 - jx1) # end-index in npick array
        if brcids[pix1]==qc1 # regular: event 1 in cluster 1, event 2 in cluster 2
            dqlat1[ix1:ix2] .= brlats[pix1]-btlats[qc1]
            dqlat2[ix1:ix2] .= brlats[pix2]-btlats[qc2]
            dqlon1[ix1:ix2] .= brlons[pix1]-btlons[qc1]
            dqlon2[ix1:ix2] .= brlons[pix2]-btlons[qc2]
            dqdep1[ix1:ix2] .= brdeps[pix1]-btdeps[qc1]
            dqdep2[ix1:ix2] .= brdeps[pix2]-btdeps[qc2]
            dqorg1[ix1:ix2] .= brorgs[pix1]-btorgs[qc1]
            dqorg2[ix1:ix2] .= brorgs[pix2]-btorgs[qc2]
            tdif21[ix1:ix2] .= tdif[jx1:jx2] # keep the tdif with same sign
        else # flipped: event 1 in cluster 2, event 2 in cluster 1
            dqlat1[ix1:ix2] .= brlats[pix2]-btlats[qc1]
            dqlat2[ix1:ix2] .= brlats[pix1]-btlats[qc2]
            dqlon1[ix1:ix2] .= brlons[pix2]-btlons[qc1]
            dqlon2[ix1:ix2] .= brlons[pix1]-btlons[qc2]
            dqdep1[ix1:ix2] .= brdeps[pix2]-btdeps[qc1]
            dqdep2[ix1:ix2] .= brdeps[pix1]-btdeps[qc2]
            dqorg1[ix1:ix2] .= brorgs[pix2]-btorgs[qc1]
            dqorg2[ix1:ix2] .= brorgs[pix1]-btorgs[qc2]
            tdif21[ix1:ix2] .= -tdif[jx1:jx2] # flip the tdif
        end
        slat21[ix1:ix2] .= slat[jx1:jx2] # stays same
        slon21[ix1:ix2] .= slon[jx1:jx2] # stays same
        phase21[ix1:ix2] .= iphase[jx1:jx2] # stays same
        ix1 = ix2 + 1 # update start index in npick array
    end

    # unweighted cluster centroid
    clat0 = (btlats[qc1] + btlats[qc2]) / 2.0
    clon0 = (btlons[qc1] + btlons[qc2]) / 2.0
    cdep0 = (btdeps[qc1] + btdeps[qc2]) / 2.0

    # run difclust (relocation norms 1, 2, 3)
    if irelonorm == 1
        clat1, clon1, cdep1, clat2, clon2, cdep2, 
        cdist, torgdif, resid, rms, rmed, resol = difclust1(
        clat0,clon0,cdep0,tdif21,phase21, slat21, slon21,
        dqlat1, dqlon1, dqdep1, dqorg1, dqlat2, dqlon2, dqdep2, dqorg2,
        ttTABs,boxwid,nit,degkm)
    elseif irelonorm == 2
        clat1, clon1, cdep1, clat2, clon2, cdep2, 
        cdist, torgdif, resid, rms, rmed, resol = difclust2(
        clat0,clon0,cdep0,tdif21,phase21, slat21, slon21,
        dqlat1, dqlon1, dqdep1, dqorg1, dqlat2, dqlon2, dqdep2, dqorg2,
        ttTABs,boxwid,nit,degkm)
    else
        clat1, clon1, cdep1, clat2, clon2, cdep2, 
        cdist, torgdif, resid, rms, rmed, resol = difclust3(
            clat0,clon0,cdep0,tdif21,phase21, slat21, slon21,
            dqlat1, dqlon1, dqdep1, dqorg1, dqlat2, dqlon2, dqdep2, dqorg2,
            ttTABs,boxwid,nit,degkm)
    end

    # careful with cluster mergers near surface
    if (min(cdep1,cdep2) < cdepmin)
        continue
    end

    # reject cluster merger if rms or median absolute residual too large or relocated dist to far
    if ((rms > rmsmax) | (rmed > rmedmax) | (cdist > distmax2))
        continue
    end

    # for robustness
    if abs(torgdif > torgdifmax)
        println("LARGE ORIGIN TIME CORRECTION: $torgdif")
        println("Possible xcor data or event list error!")
        continue
        #exit()
    end

    # fraction of events in each cluster
    fracC1 = Float64(nb1)/Float64(nb1+nb2)
    fracC2 = 1.0-fracC1

    # original centroid of combined cluster
    #creflat = cosd(clat0) # needed to scale dlon
    cxlat00 = btlats[qc1]*fracC1+btlats[qc2]*fracC2
    cxlon00 = btlons[qc1]*fracC1+btlons[qc2]*fracC2
    cxdep00 = btdeps[qc1]*fracC1+btdeps[qc2]*fracC2

    # new centroid of combined cluster
    cxlat11 = clat1*fracC1+clat2*fracC2
    cxlon11 = clon1*fracC1+clon2*fracC2
    cxdep11 = cdep1*fracC1+cdep2*fracC2

    # offset between two (b/c not all links used)
    dcxlat = cxlat11-cxlat00
    dcxlon = cxlon11-cxlon00
    dcxdep = cxdep11-cxdep00


    # check relative shift of cluster 1 (subtracting possible DC offset)
    # offsets for cluster 1: new loc - old loc
    qlat_off1 = clat1 - btlats[qc1]
    qlon_off1 = clon1 - btlons[qc1]
    qdep_off1 = cdep1 - btdeps[qc1]        
    if abs(qdep_off1-dcxdep) > vshiftmax
        continue
    elseif ((qlon_off1-dcxlon)*cosd(clat0))^2 + (qlat_off1-dcxlat)^2 > hshiftmaxD2 # in squared degrees
        continue
    end 

    # check relative shift of cluster 2 (subtracting possible DC offset)
    # offsets for cluster 2: new loc - old loc
    qlat_off2 = clat2 - btlats[qc2]
    qlon_off2 = clon2 - btlons[qc2]
    qdep_off2 = cdep2 - btdeps[qc2]
    if abs(qdep_off2 - dcxdep) > vshiftmax
        continue
    elseif ((qlon_off2-dcxlon)*cosd(clat0))^2 + (qlat_off2-dcxlat)^2 > hshiftmaxD2 # in squared degrees
        continue
    end 

    # ok, cluster merger is approved!
    #println("Approved!")

    # origin time updates
    @inbounds begin
    #qtim_off1 = btorgs[qc1] - torgdif/2.0 # symmetric shift "backward" if positive
    #qtim_off2 = btorgs[qc2] + torgdif/2.0 # symmetric shift "forward" if positive
    qtim_off1 = -torgdif/2.0 # symmetric shift "backward" if positive (torg always 0)
    qtim_off2 =  torgdif/2.0 # symmetric shift "forward" if positive (torg always 0)
    cxorg11 = fracC1*qtim_off1 + fracC2*qtim_off2

    # update locations: cluster 1
    brlats[c1qixs] .+= qlat_off1
    brlons[c1qixs] .+= qlon_off1
    brdeps[c1qixs] .+= qdep_off1
    brorgs[c1qixs] .+= qtim_off1

    # update locations: cluster 2
    brlats[c2qixs] .+= qlat_off2
    brlons[c2qixs] .+= qlon_off2
    brdeps[c2qixs] .+= qdep_off2
    brorgs[c2qixs] .+= qtim_off2

    # evacuate tree 2, assign events to tree 1
    brcids[c2qixs] .= qc1
    btnbranch[qc1] += nb2
    #btnbranch[qc2] = 0  # correct, but not needed

    end # end of @inbounds

    # merged cluster
    union!(cid2pairD[qc1],cid2pairD[qc2])
    union!(cid2qixD[qc1],c2qixs)


    # # align cluster and catalog centroids, set average time to zero
    # #  (need to do this b/c not all events used in difclust
    iclust = cid2qixD[qc1]
    @inbounds begin
    brorgs[iclust] .-= cxorg11 #mean(brorgs[iclust])
    brlats[iclust] .+= (cxlat00 - cxlat11)
    btlats[iclust] .= cxlat00
    brlons[iclust] .+= (cxlon00 - cxlon11)
    btlons[iclust] .= cxlon00
    brdeps[iclust] .+= (cxdep00 - cxdep11)
    btdeps[iclust] .= cxdep00
    #btorgs[iclust] .= 0.0 # always zero, never updated
    end

    end # ---- end of clustering loop

    # return results
    bnb = btnbranch[brcids]
    return brlats, brlons, brdeps, brorgs, brcids, bnb
end

# function with i64 ixx arrays
function clustertree(pqix1::Vector{Int32},pqix2::Vector{Int32},ixx1::Vector{Int64},ixx2::Vector{Int64},
    tdif::Vector{Float64},slat::Vector{Float64},slon::Vector{Float64},iphase::Vector{Int8},
    qlats::Vector{Float64},qlons::Vector{Float64},qdeps::Vector{Float64},ttTABs,
    nit::Int64,boxwid::Float64,degkm::Float64,irelonorm::Int64,rmsmax::Float64,rmedmax::Float64,
    distmax::Float64,distmax2::Float64,hshiftmax::Float64,vshiftmax::Float64,torgdifmax::Float64,
    nupdate::Int64,maxlink::Int64)

    # setup parameters
    cdepmin = min(0.0,minimum(qdeps))
    hshiftmaxD2 = (hshiftmax/degkm)^2 # in squared degrees
    distmax22 = distmax^2 # in squared km

    # array sizes
    nq = length(qlats)
    bnpair = length(pqix1)

    # initialize relocated event arrays
    brlats = copy(qlats)
    brlons = copy(qlons)
    brdeps = copy(qdeps)
    brorgs = zeros(nq)
    brcids = Vector{Int32}(1:nq)

    # initialize clustering tree arrays
    btlats = copy(brlats)
    btlons = copy(brlons)
    btdeps = copy(brdeps)
    btorgs = zeros(nq)
    btnbranch = ones(Int32,nq)

    # dictionary with cid => event indices (note cid=qix on initialization)
    cid2qixD = Dict(qix=>[qix] for qix in brcids)

    # dictionary with cid => event pair indices
    cid2pairD = Dict( qix => Vector{Int32}(findall(
    (pqix1.==qix).|(pqix2 .== qix))) for qix in brcids) # cid=qix on initialization


    # loop over event pairs
    @inbounds for ip in Vector{Int32}(1:bnpair)

    # Progress
    if (mod(ip,nupdate)==0)
        println("Thread: ", Threads.threadid(), "--> Working on sorted pair: $ip/$bnpair")
    end

    # get event pair information
    qix1, qix2 = pqix1[ip], pqix2[ip]
    qc1, qc2 = brcids[qix1], brcids[qix2]
    if qc1 == qc2; continue; end # skip if in same cluster
    nb1, nb2 = btnbranch[qc1], btnbranch[qc2]

    # check to see if clusters are too far apart
    cdist22 = (btdeps[qc2]-btdeps[qc1])^2 + map_distance(
    btlats[qc1],btlons[qc1],btlats[qc2],btlons[qc2])^2
    if cdist22 > distmax22
        continue
    end

    # find event ids in each cluster
    c1qixs = cid2qixD[qc1] 
    c2qixs = cid2qixD[qc2]

    # only one link
    if (nb1==1)&(nb2==1)

        # only one link
        nlink = 1
        linx = [ip]

    # find all links between clusters
    else

        # # event pairs linking clusters
        # #   - pairs above are either processed and in same cluster
        # #     or caused issues when trying to link
        cpix1, cpix2 = cid2pairD[qc1],cid2pairD[qc2]
        @views linx = intersect(cpix1[cpix1.>=ip],cpix2[cpix2.>=ip])

        # keeping only best linx
        nlink = length(linx)
        if nlink > maxlink
            linx = partialsort(linx,1:maxlink)
        end 

    end

    # count number of picks
    npick = sum([ixx2[jj] - ixx1[jj] + 1 for jj in linx])

    # extract locations relative to centroid
    dqlat1, dqlon1 = zeros(npick),zeros(npick)
    dqdep1, dqorg1 = zeros(npick),zeros(npick)
    dqlat2, dqlon2 = zeros(npick),zeros(npick)
    dqdep2, dqorg2 = zeros(npick),zeros(npick)
    phase21, slat21, slon21 = zeros(Int8,npick), zeros(npick), zeros(npick)
    tdif21 = zeros(npick)
    ix1 = 1 # start index for event pair in the arrays above
    @inbounds for jj in linx
        pix1, pix2, jx1, jx2 = pqix1[jj],pqix2[jj],ixx1[jj],ixx2[jj]
        ix2 = ix1 + (jx2 - jx1) # end-index in npick array
        if brcids[pix1]==qc1 # regular: event 1 in cluster 1, event 2 in cluster 2
            dqlat1[ix1:ix2] .= brlats[pix1]-btlats[qc1]
            dqlat2[ix1:ix2] .= brlats[pix2]-btlats[qc2]
            dqlon1[ix1:ix2] .= brlons[pix1]-btlons[qc1]
            dqlon2[ix1:ix2] .= brlons[pix2]-btlons[qc2]
            dqdep1[ix1:ix2] .= brdeps[pix1]-btdeps[qc1]
            dqdep2[ix1:ix2] .= brdeps[pix2]-btdeps[qc2]
            dqorg1[ix1:ix2] .= brorgs[pix1]-btorgs[qc1]
            dqorg2[ix1:ix2] .= brorgs[pix2]-btorgs[qc2]
            tdif21[ix1:ix2] .= tdif[jx1:jx2] # keep the tdif with same sign
        else # flipped: event 1 in cluster 2, event 2 in cluster 1
            dqlat1[ix1:ix2] .= brlats[pix2]-btlats[qc1]
            dqlat2[ix1:ix2] .= brlats[pix1]-btlats[qc2]
            dqlon1[ix1:ix2] .= brlons[pix2]-btlons[qc1]
            dqlon2[ix1:ix2] .= brlons[pix1]-btlons[qc2]
            dqdep1[ix1:ix2] .= brdeps[pix2]-btdeps[qc1]
            dqdep2[ix1:ix2] .= brdeps[pix1]-btdeps[qc2]
            dqorg1[ix1:ix2] .= brorgs[pix2]-btorgs[qc1]
            dqorg2[ix1:ix2] .= brorgs[pix1]-btorgs[qc2]
            tdif21[ix1:ix2] .= -tdif[jx1:jx2] # flip the tdif
        end
        slat21[ix1:ix2] .= slat[jx1:jx2] # stays same
        slon21[ix1:ix2] .= slon[jx1:jx2] # stays same
        phase21[ix1:ix2] .= iphase[jx1:jx2] # stays same
        ix1 = ix2 + 1 # update start index in npick array
    end

    # unweighted cluster centroid
    clat0 = (btlats[qc1] + btlats[qc2]) / 2.0
    clon0 = (btlons[qc1] + btlons[qc2]) / 2.0
    cdep0 = (btdeps[qc1] + btdeps[qc2]) / 2.0

    # run difclust (relocation norms 1, 2, 3)
    if irelonorm == 1
        clat1, clon1, cdep1, clat2, clon2, cdep2, 
        cdist, torgdif, resid, rms, rmed, resol = difclust1(
        clat0,clon0,cdep0,tdif21,phase21, slat21, slon21,
        dqlat1, dqlon1, dqdep1, dqorg1, dqlat2, dqlon2, dqdep2, dqorg2,
        ttTABs,boxwid,nit,degkm)
    elseif irelonorm == 2
        clat1, clon1, cdep1, clat2, clon2, cdep2, 
        cdist, torgdif, resid, rms, rmed, resol = difclust2(
        clat0,clon0,cdep0,tdif21,phase21, slat21, slon21,
        dqlat1, dqlon1, dqdep1, dqorg1, dqlat2, dqlon2, dqdep2, dqorg2,
        ttTABs,boxwid,nit,degkm)
    else
        clat1, clon1, cdep1, clat2, clon2, cdep2, 
        cdist, torgdif, resid, rms, rmed, resol = difclust3(
            clat0,clon0,cdep0,tdif21,phase21, slat21, slon21,
            dqlat1, dqlon1, dqdep1, dqorg1, dqlat2, dqlon2, dqdep2, dqorg2,
            ttTABs,boxwid,nit,degkm)
    end

    # careful with cluster mergers near surface
    if (min(cdep1,cdep2) < cdepmin)
        continue
    end

    # for robustness
    if abs(torgdif > torgdifmax)
        println("LARGE ORIGIN TIME CORRECTION: $torgdif")
        println("Likely xcor data or event list error!")
        exit()
    end

    # reject cluster merger if rms or median absolute residual too large or relocated dist to far
    if ((rms > rmsmax) | (rmed > rmedmax) | (cdist > distmax2))
        continue
    end

    # fraction of events in each cluster
    fracC1 = Float64(nb1)/Float64(nb1+nb2)
    fracC2 = 1.0-fracC1

    # original centroid of combined cluster
    #creflat = cosd(clat0) # needed to scale dlon
    cxlat00 = btlats[qc1]*fracC1+btlats[qc2]*fracC2
    cxlon00 = btlons[qc1]*fracC1+btlons[qc2]*fracC2
    cxdep00 = btdeps[qc1]*fracC1+btdeps[qc2]*fracC2

    # new centroid of combined cluster
    cxlat11 = clat1*fracC1+clat2*fracC2
    cxlon11 = clon1*fracC1+clon2*fracC2
    cxdep11 = cdep1*fracC1+cdep2*fracC2

    # offset between two (b/c not all links used)
    dcxlat = cxlat11-cxlat00
    dcxlon = cxlon11-cxlon00
    dcxdep = cxdep11-cxdep00


    # check relative shift of cluster 1 (subtracting possible DC offset)
    # offsets for cluster 1: new loc - old loc
    qlat_off1 = clat1 - btlats[qc1]
    qlon_off1 = clon1 - btlons[qc1]
    qdep_off1 = cdep1 - btdeps[qc1]        
    if abs(qdep_off1-dcxdep) > vshiftmax
        continue
    elseif ((qlon_off1-dcxlon)*cosd(clat0))^2 + (qlat_off1-dcxlat)^2 > hshiftmaxD2 # in squared degrees
        continue
    end 

    # check relative shift of cluster 2 (subtracting possible DC offset)
    # offsets for cluster 2: new loc - old loc
    qlat_off2 = clat2 - btlats[qc2]
    qlon_off2 = clon2 - btlons[qc2]
    qdep_off2 = cdep2 - btdeps[qc2]
    if abs(qdep_off2 - dcxdep) > vshiftmax
        continue
    elseif ((qlon_off2-dcxlon)*cosd(clat0))^2 + (qlat_off2-dcxlat)^2 > hshiftmaxD2 # in squared degrees
        continue
    end 

    # ok, cluster merger is approved!
    #println("Approved!")

    # origin time updates
    @inbounds begin
    #qtim_off1 = btorgs[qc1] - torgdif/2.0 # symmetric shift "backward" if positive
    #qtim_off2 = btorgs[qc2] + torgdif/2.0 # symmetric shift "forward" if positive
    qtim_off1 = -torgdif/2.0 # symmetric shift "backward" if positive (torg always 0)
    qtim_off2 =  torgdif/2.0 # symmetric shift "forward" if positive (torg always 0)
    cxorg11 = fracC1*qtim_off1 + fracC2*qtim_off2

    # update locations: cluster 1
    brlats[c1qixs] .+= qlat_off1
    brlons[c1qixs] .+= qlon_off1
    brdeps[c1qixs] .+= qdep_off1
    brorgs[c1qixs] .+= qtim_off1

    # update locations: cluster 2
    brlats[c2qixs] .+= qlat_off2
    brlons[c2qixs] .+= qlon_off2
    brdeps[c2qixs] .+= qdep_off2
    brorgs[c2qixs] .+= qtim_off2

    # evacuate tree 2, assign events to tree 1
    brcids[c2qixs] .= qc1
    btnbranch[qc1] += nb2
    #btnbranch[qc2] = 0  # correct, but not needed

    end # end of @inbounds

    # merged cluster
    union!(cid2pairD[qc1],cid2pairD[qc2])
    union!(cid2qixD[qc1],c2qixs)


    # # align cluster and catalog centroids, set average time to zero
    # #  (need to do this b/c not all events used in difclust
    iclust = cid2qixD[qc1]
    @inbounds begin
    brorgs[iclust] .-= cxorg11 #mean(brorgs[iclust])
    brlats[iclust] .+= (cxlat00 - cxlat11)
    btlats[iclust] .= cxlat00
    brlons[iclust] .+= (cxlon00 - cxlon11)
    btlons[iclust] .= cxlon00
    brdeps[iclust] .+= (cxdep00 - cxdep11)
    btdeps[iclust] .= cxdep00
    #btorgs[iclust] .= 0.0 # always zero, never updated
    end

    end # ---- end of clustering loop

    # return results
    bnb = btnbranch[brcids]
    return brlats, brlons, brdeps, brorgs, brcids, bnb
end