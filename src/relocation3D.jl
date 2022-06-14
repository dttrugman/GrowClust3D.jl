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

# Inputs: qX0  =  reference center point xpos (km)
#         qY0  =  reference center point ypos (km)
#         qZ0  =  reference center point depth (km)
#         tdif =  array (len=npick) of dif times, t2-t1 (s)
#         itab =  array (len=npick) with table index numbers
#         qX1  =  array (len=npick) of events in cluster1 xpos offsets from centroid
#         qY1  =  array (len=npick) of events in cluster1 ypos offsets
#         qZ1  =  array (len=npick) of events in cluster1 depth offsets
#         qT1  =  array (len=npick) of events in cluster1 time offsets
#         qX2  =  array (len=npick) of events in cluster2 xpos offsets from centroid
#         qY2  =  array (len=npick) of events in cluster2 ypos offsets
#         qZ2  =  array (len=npick) of events in cluster2 depth offsets
#         qT2  =  array (len=npick) of events in cluster2 time offsets
#         ttTABs =  travel time tables for different phases and stations
#         boxwid =  starting box width (km)
#         nit    =  number of iterations to perform
# Returns: cX1  =  best-fitting xpos (km) for first cluster centroid
#          cY1  =  best-fitting ypos (km) for first cluster centroid
#          cZ1  =  best-fitting depth (km) of first cluster centroid
#          cX2  =  best-fitting xpos (km) for second cluster centroid
#          cY2  =  best-fitting ypos (km) for second cluster centroid
#          cZ2  =  best-fitting depth (km) of second cluster centroid
#          cdist  =  cluster separation distance (km)
#          torgdif=  origin time difference, i.e., t2-t1 median residual (>0 when 2 is later than 1)
#          resid  =  array (len=npick) of residuals (s) between observed tdif (tt) and predicted
#          rms    =  rms residual ( sqrt( sum(resid**2)/npick) )
#          rmed   =  median absolute value residual ( median( abs(resid) ) )
#          resol  =  nominal resolution (m) of final box

##### difclust1: L1 residual norm
function difclust1_3D(qX0::Float64,qY0::Float64,qZ0::Float64,
    tdif::Vector{Float32},itab::Vector{Int16},
    qX1::Vector{Float64},qY1::Vector{Float64},
    qZ1::Vector{Float64},qT1::Vector{Float32},
    qX2::Vector{Float64},qY2::Vector{Float64},
    qZ2::Vector{Float64},qT2::Vector{Float32},
    ttTABs::AbstractVector{},boxwid::Float64,nit::Int64)

# initialize variables to track best solution
fxbest1, fybest1, fzbest1 = 0.0, 0.0, 0.0
fxbest2, fybest2, fzbest2 = 0.0, 0.0, 0.0
torgdif = Float32(0.0)

# extract npick
npick = length(tdif)
resid = zeros(Float32,npick)

# initialize box
dX0, dY0, dZ0 = 0.0, 0.0, 0.0 # C1 shift from initial centroid
dX, dY = 0.5*boxwid, 0.5*boxwid # current box bounds

# careful near surface (added 07/2018)
zboxwid = boxwid
if (qZ0 < 0.0)
    zboxwid = min(1.0, boxwid)
elseif (qZ0 < boxwid/2.0)
    zboxwid = max(1.0,2.0*qZ0)
end
dZ = 0.5*zboxwid # analagous to dX and dY

# start iterations
for it=1:nit
    
    # grid search over box for best fit
    fitbest = Float32(1.0e20)
    tbest = Float32(0.0)
    
    # get trial locations: f1, f2 (3x3 grid for this box iteration)
    for iy = -1.0:1.0
        fY1 = qY0 + dY0 + dY*iy # centroid + shift + trial offset
        fY2 = qY0 - dY0 - dY*iy # centroid - shift - trial offset
        for ix = -1.0:1.0
            fX1 = qX0 + dX0 + dX*ix  # centroid + shift + trial offset
            fX2 = qX0 - dX0 - dX*ix # centroid - shift - trial offset
            for iz = -1.0:1.0
                fZ1 = qZ0 + dZ0 + dZ*iz  # centroid + shift + trial offset
                fZ2 = qZ0 - dZ0 - dZ*iz # centroid - shift - trial offset
                
                # compute predicted travel time and residuals w/observed           
                @inbounds for ii=1:npick
                    tt1 = Float32(ttTABs[itab[ii]](qX1[ii]+fX1,qY1[ii]+fY1,qZ1[ii]+fZ1))
                    tt2 = Float32(ttTABs[itab[ii]](qX2[ii]+fX2,qY2[ii]+fY2,qZ2[ii]+fZ2))
                    pdif = (tt2 + qT2[ii]) - (tt1 + qT1[ii])
                    resid[ii] = tdif[ii] - pdif # accounts for otime adjustment
                end
                
                # compute fit: L1 norm
                residval = median(resid)
                fit = sum(abs.(resid.-residval))
                
                # update best fit
                if fit < fitbest
                    fitbest = fit
                    fxbest1 = fX1
                    fybest1 = fY1
                    fzbest1 = fZ1
                    fxbest2 = fX2
                    fybest2 = fY2
                    fzbest2 = fZ2
                    tbest = residval
                end
                
            end
        end
    end # end of grid search
    
    # update best shift
    dX0 = fxbest1 - qX0
    dY0 = fybest1 - qY0
    dZ0 = fzbest1 - qZ0
    torgdif = tbest # otime difference: residual between observed and predicted tdif

    # shrink box by 2/3 each iteration
    dX *= 2.0/3.0      
    dY *= 2.0/3.0
    dZ *= 2.0/3.0
    
end # end loop over iterations

# output final locations
cX1 = fxbest1
cY1 = fybest1
cZ1 = fzbest1
cX2 = fxbest2
cY2 = fybest2
cZ2 = fzbest2
resol = dY/(2.0/3.0)

# compute distance between cluster centroids
cdist = sqrt((cX2-cX1)^2 + (cY2-cY1)^2 + (cZ2-cZ1)^2)
                
# compute residual between observed and predicted travel time
@inbounds for ii=1:npick
    tt1 = Float32(ttTABs[itab[ii]](cX1+qX1[ii],cY1+qY1[ii],cZ1+qZ1[ii]))
    tt2 = Float32(ttTABs[itab[ii]](cX2+qX2[ii],cY2+qY2[ii],cZ2+qZ2[ii]))
    pdif = (tt2 + qT2[ii]) - (tt1 + qT1[ii])
    resid[ii] = tdif[ii] - pdif - torgdif
end

# rms residual and median absolute residual
rms = sqrt.(mean(resid.^2))
rmed = mad(resid,normalize=false) # default is to equate to sdev

# return results
return cX1, cY1, cZ1, cX2, cY2, cZ2, cdist, torgdif, resid, rms, rmed, resol

end

##### difclust2: L2 residual norm
function difclust2_3D(qX0::Float64,qY0::Float64,qZ0::Float64,
    tdif::Vector{Float32},itab::Vector{Int16},
    qX1::Vector{Float64},qY1::Vector{Float64},
    qZ1::Vector{Float64},qT1::Vector{Float32},
    qX2::Vector{Float64},qY2::Vector{Float64},
    qZ2::Vector{Float64},qT2::Vector{Float32},
    ttTABs::AbstractVector{},boxwid::Float64,nit::Int64)

# initialize variables to track best solution
fxbest1, fybest1, fzbest1 = 0.0, 0.0, 0.0
fxbest2, fybest2, fzbest2 = 0.0, 0.0, 0.0
torgdif = Float32(0.0)

# extract npick
npick = length(tdif)
resid = zeros(Float32,npick)

# initialize box
dX0, dY0, dZ0 = 0.0, 0.0, 0.0 # C1 shift from initial centroid
dX, dY = 0.5*boxwid, 0.5*boxwid # current box bounds

# careful near surface (added 07/2018)
zboxwid = boxwid
if (qZ0 < 0.0)
    zboxwid = min(1.0, boxwid)
elseif (qZ0 < boxwid/2.0)
    zboxwid = max(1.0,2.0*qZ0)
end
dZ = 0.5*zboxwid # analagous to dX and dY

# start iterations
for it=1:nit
    
    # grid search over box for best fit
    fitbest = Float32(1.0e20)
    tbest = Float32(0.0)
    
    # get trial locations: f1, f2 (3x3 grid for this box iteration)
    for iy = -1.0:1.0
        fY1 = qY0 + dY0 + dY*iy # centroid + shift + trial offset
        fY2 = qY0 - dY0 - dY*iy # centroid - shift - trial offset
        for ix = -1.0:1.0
            fX1 = qX0 + dX0 + dX*ix  # centroid + shift + trial offset
            fX2 = qX0 - dX0 - dX*ix # centroid - shift - trial offset
            for iz = -1.0:1.0
                fZ1 = qZ0 + dZ0 + dZ*iz  # centroid + shift + trial offset
                fZ2 = qZ0 - dZ0 - dZ*iz # centroid - shift - trial offset
                
                # compute predicted travel time and residuals w/observed             
                @inbounds for ii=1:npick
                    tt1 = Float32(ttTABs[itab[ii]](qX1[ii]+fX1,qY1[ii]+fY1,qZ1[ii]+fZ1))
                    tt2 = Float32(ttTABs[itab[ii]](qX2[ii]+fX2,qY2[ii]+fY2,qZ2[ii]+fZ2))
                    pdif = (tt2 + qT2[ii]) - (tt1 + qT1[ii])
                    resid[ii] = tdif[ii] - pdif # accounts for otime adjustment
                end
                
                # compute fit: L2 norm
                residval = mean(resid)
                fit = sum((resid.-residval).^2)
                
                # update best fit
                if fit < fitbest
                    fitbest = fit
                    fxbest1 = fX1
                    fybest1 = fY1
                    fzbest1 = fZ1
                    fxbest2 = fX2
                    fybest2 = fY2
                    fzbest2 = fZ2
                    tbest = residval
                end
                
            end
        end
    end # end of grid search
    
    # update best shift
    dX0 = fxbest1 - qX0
    dY0 = fybest1 - qY0
    dZ0 = fzbest1 - qZ0
    torgdif = tbest # otime difference: residual between observed and predicted tdif

    # shrink box by 2/3 each iteration
    dX *= 2.0/3.0      
    dY *= 2.0/3.0
    dZ *= 2.0/3.0
    
end # end loop over iterations

# output final locations
cX1 = fxbest1
cY1 = fybest1
cZ1 = fzbest1
cX2 = fxbest2
cY2 = fybest2
cZ2 = fzbest2
resol = dY/(2.0/3.0)

# compute distance between cluster centroids
cdist = sqrt((cX2-cX1)^2 + (cY2-cY1)^2 + (cZ2-cZ1)^2)
                
# compute residual between observed and predicted travel time
@inbounds for ii=1:npick
    tt1 = Float32(ttTABs[itab[ii]](cX1+qX1[ii],cY1+qY1[ii],cZ1+qZ1[ii]))
    tt2 = Float32(ttTABs[itab[ii]](cX2+qX2[ii],cY2+qY2[ii],cZ2+qZ2[ii]))
    pdif = (tt2 + qT2[ii]) - (tt1 + qT1[ii])
    resid[ii] = tdif[ii] - pdif - torgdif
end

# rms residual and median absolute residual
rms = sqrt.(mean(resid.^2))
rmed = mad(resid,normalize=false) # default is to equate to sdev

# return results
return cX1, cY1, cZ1, cX2, cY2, cZ2, cdist, torgdif, resid, rms, rmed, resol

end

##### difclust3: L3 / Robomean Norm
function difclust3_3D(qX0::Float64,qY0::Float64,qZ0::Float64,
    tdif::Vector{Float32},itab::Vector{Int16}, 
    qX1::Vector{Float64},qY1::Vector{Float64},
    qZ1::Vector{Float64},qT1::Vector{Float32},
    qX2::Vector{Float64},qY2::Vector{Float64},
    qZ2::Vector{Float64},qT2::Vector{Float32},
    ttTABs::AbstractVector{},boxwid::Float64,nit::Int64)

# initialize variables to track best solution
fxbest1, fybest1, fzbest1 = 0.0, 0.0, 0.0
fxbest2, fybest2, fzbest2 = 0.0, 0.0, 0.0
torgdif = Float32(0.0)

# extract npick
npick = length(tdif)
resid = zeros(Float32,npick)

# initialize box
dX0, dY0, dZ0 = 0.0, 0.0, 0.0 # C1 shift from initial centroid
dX, dY = 0.5*boxwid, 0.5*boxwid # current box bounds

# careful near surface (added 07/2018)
zboxwid = boxwid
if (qZ0 < 0.0)
    zboxwid = min(1.0, boxwid)
elseif (qZ0 < boxwid/2.0)
    zboxwid = max(1.0,2.0*qZ0)
end
dZ = 0.5*zboxwid # analagous to dX and dY

# start iterations
for it=1:nit
    
    # grid search over box for best fit
    fitbest = Float32(1.0e20)
    tbest = Float32(0.0)
    
    # get trial locations: f1, f2 (3x3 grid for this box iteration)
    for iy = -1.0:1.0
        fY1 = qY0 + dY0 + dY*iy # centroid + shift + trial offset
        fY2 = qY0 - dY0 - dY*iy # centroid - shift - trial offset
        for ix = -1.0:1.0
            fX1 = qX0 + dX0 + dX*ix  # centroid + shift + trial offset
            fX2 = qX0 - dX0 - dX*ix # centroid - shift - trial offset
            for iz = -1.0:1.0
                fZ1 = qZ0 + dZ0 + dZ*iz  # centroid + shift + trial offset
                fZ2 = qZ0 - dZ0 - dZ*iz # centroid - shift - trial offset
                
                # compute predicted travel time and residuals w/observed          
                @inbounds for ii=1:npick
                    tt1 = Float32(ttTABs[itab[ii]](qX1[ii]+fX1,qY1[ii]+fY1,qZ1[ii]+fZ1))
                    tt2 = Float32(ttTABs[itab[ii]](qX2[ii]+fX2,qY2[ii]+fY2,qZ2[ii]+fZ2))
                    pdif = (tt2 + qT2[ii]) - (tt1 + qT1[ii])
                    resid[ii] = tdif[ii] - pdif # accounts for otime adjustment
                end
                
                # compute fit: L3 norm
                residval, fit = robomean(resid,0.1,10)
                
                # update best fit
                if fit < fitbest
                    fitbest = fit
                    fxbest1 = fX1
                    fybest1 = fY1
                    fzbest1 = fZ1
                    fxbest2 = fX2
                    fybest2 = fY2
                    fzbest2 = fZ2
                    tbest = residval
                end
                
            end
        end
    end # end of grid search
    
    # update best shift
    dX0 = fxbest1 - qX0
    dY0 = fybest1 - qY0
    dZ0 = fzbest1 - qZ0
    torgdif = tbest # otime difference: residual between observed and predicted tdif

    # shrink box by 2/3 each iteration
    dX *= 2.0/3.0      
    dY *= 2.0/3.0
    dZ *= 2.0/3.0
    
end # end loop over iterations

# output final locations
cX1 = fxbest1
cY1 = fybest1
cZ1 = fzbest1
cX2 = fxbest2
cY2 = fybest2
cZ2 = fzbest2
resol = dY/(2.0/3.0)

# compute distance between cluster centroids
cdist = sqrt((cX2-cX1)^2 + (cY2-cY1)^2 + (cZ2-cZ1)^2)
                
# compute residual between observed and predicted travel time
@inbounds for ii=1:npick
    tt1 = Float32(ttTABs[itab[ii]](cX1+qX1[ii],cY1+qY1[ii],cZ1+qZ1[ii]))
    tt2 = Float32(ttTABs[itab[ii]](cX2+qX2[ii],cY2+qY2[ii],cZ2+qZ2[ii]))
    pdif = (tt2 + qT2[ii]) - (tt1 + qT1[ii])
    resid[ii] = tdif[ii] - pdif - torgdif
end

# rms residual and median absolute residual
rms = sqrt.(mean(resid.^2))
rmed = mad(resid,normalize=false) # default is to equate to sdev

# return results
return cX1, cY1, cZ1, cX2, cY2, cZ2, cdist, torgdif, resid, rms, rmed, resol

end

#####################################

### CLUSTERTREE: Runs the the main growclust algorithm, relocated events using xcorr data / clustering.
#
# Inputs:
# - pqix1, pqix2: vectors of serial event numbers for each event pair
# - ixx1, ixx2: start and stop indices in xcorr arrays for each pair
# - tdif, itab: differential travel times and table index for each xcorr data
# - qXs, qXs, qZs: initial event locations
# - ttTABS: array of travel time interpolation objects
# - nit, boxwid, degkm, irelonorm: difclust parameters
# - rmsmax, rmedmax, distmax, distmax2, hshiftmax, vshifmax, torgdifmax: cluster merge parameters
# - nupdate: counter for iteration update
# - maxlink: maximum number of event pairs used to link clusters
#
# Returns:
# - brXs, brYs, brZs: relocated event locations
# - brcids, bnb: cluster ids and number of events in each cluster 
#
# function with i32 ixx arrays
function clustertree_3D(pqix1::Vector{Int32},pqix2::Vector{Int32},ixx1::Vector{Int32},ixx2::Vector{Int32},
    tdif::Vector{Float32},itab::Vector{Int16},
    qXs::Vector{Float64},qYs::Vector{Float64},qZs::Vector{Float64},ttTABs::AbstractVector{},
    nit::Int64,boxwid::Float64,irelonorm::Int64,rmsmax::Float32,rmedmax::Float32,
    distmax::Float64,distmax2::Float64,hshiftmax::Float64,vshiftmax::Float64,torgdifmax::Float32,
    nupdate::Int64,maxlink::Int64)

    # setup parameters
    cZmin = min(0.0,minimum(qZs))
    distmax22 = distmax^2 # in squared km

    # array sizes
    nq = length(qXs)
    bnpair = length(pqix1)

    # initialize relocated event arrays
    brXs = copy(qXs)
    brYs = copy(qYs)
    brZs = copy(qZs)
    brorgs = zeros(Float32,nq)
    brcids = Vector{Int32}(1:nq)

    # initialize clustering tree arrays
    btXs = copy(brXs)
    btYs = copy(brYs)
    btZs = copy(brZs)
    btorgs = zeros(Float32,nq)
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
    cdist22 = (btXs[qc2]-btXs[qc1])^2 + (btYs[qc2]-btYs[qc1])^2 + (btZs[qc2]-btZs[qc1])^2
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
    dqX1, dqY1 = zeros(npick),zeros(npick)
    dqZ1, dqorg1 = zeros(npick),zeros(Float32,npick)
    dqX2, dqY2 = zeros(npick),zeros(npick)
    dqZ2, dqorg2 = zeros(npick),zeros(Float32,npick)
    tab21 = zeros(Int16,npick)
    tdif21 = zeros(Float32,npick)
    ix1 = 1 # start index for event pair in the arrays above
    @inbounds for jj in linx
        pix1, pix2, jx1, jx2 = pqix1[jj],pqix2[jj],ixx1[jj],ixx2[jj]
        ix2 = ix1 + (jx2 - jx1) # end-index in npick array
        if brcids[pix1]==qc1 # regular: event 1 in cluster 1, event 2 in cluster 2
            dqX1[ix1:ix2] .= brXs[pix1]-btXs[qc1]
            dqX2[ix1:ix2] .= brXs[pix2]-btXs[qc2]
            dqY1[ix1:ix2] .= brYs[pix1]-btYs[qc1]
            dqY2[ix1:ix2] .= brYs[pix2]-btYs[qc2]
            dqZ1[ix1:ix2] .= brZs[pix1]-btZs[qc1]
            dqZ2[ix1:ix2] .= brZs[pix2]-btZs[qc2]
            dqorg1[ix1:ix2] .= brorgs[pix1]-btorgs[qc1]
            dqorg2[ix1:ix2] .= brorgs[pix2]-btorgs[qc2]
            tdif21[ix1:ix2] .= tdif[jx1:jx2] # keep the tdif with same sign
        else # flipped: event 1 in cluster 2, event 2 in cluster 1
            dqX1[ix1:ix2] .= brXs[pix2]-btXs[qc1]
            dqX2[ix1:ix2] .= brXs[pix1]-btXs[qc2]
            dqY1[ix1:ix2] .= brYs[pix2]-btYs[qc1]
            dqY2[ix1:ix2] .= brYs[pix1]-btYs[qc2]
            dqZ1[ix1:ix2] .= brZs[pix2]-btZs[qc1]
            dqZ2[ix1:ix2] .= brZs[pix1]-btZs[qc2]
            dqorg1[ix1:ix2] .= brorgs[pix2]-btorgs[qc1]
            dqorg2[ix1:ix2] .= brorgs[pix1]-btorgs[qc2]
            tdif21[ix1:ix2] .= -tdif[jx1:jx2] # flip the tdif
        end
        tab21[ix1:ix2] .= itab[jx1:jx2] # stays same
        ix1 = ix2 + 1 # update start index in npick array
    end

    # unweighted cluster centroid
    cX0 = (btXs[qc1] + btXs[qc2]) / 2.0
    cY0 = (btYs[qc1] + btYs[qc2]) / 2.0
    cZ0 = (btZs[qc1] + btZs[qc2]) / 2.0

    # run difclust (relocation norms 1, 2, 3)
    if irelonorm == 1
        cX1, cY1, cZ1, cX2, cY2, cZ2, 
        cdist, torgdif, resid, rms, rmed, resol = difclust1_3D(
        cX0,cY0,cZ0,tdif21,tab21, 
        dqX1, dqY1, dqZ1, dqorg1,
        dqX2, dqY2, dqZ2, dqorg2,
        ttTABs,boxwid,nit)
    elseif irelonorm == 2
        cX1, cY1, cZ1, cX2, cY2, cZ2, 
        cdist, torgdif, resid, rms, rmed, resol = difclust2_3D(
        cX0,cY0,cZ0,tdif21,tab21,
        dqX1, dqY1, dqZ1, dqorg1, 
        dqX2, dqY2, dqZ2, dqorg2,
        ttTABs,boxwid,nit)
    else
        cX1, cY1, cZ1, cX2, cY2, cZ2, 
        cdist, torgdif, resid, rms, rmed, resol = difclust3_3D(
        cX0,cY0,cZ0,tdif21,tab21,
        dqX1, dqY1, dqZ1, dqorg1,
        dqX2, dqY2, dqZ2, dqorg2,
        ttTABs,boxwid,nit)
    end

    # careful with cluster mergers near surface
    if (min(cZ1,cZ2) < cZmin)
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
    cxX00 = btXs[qc1]*fracC1+btXs[qc2]*fracC2
    cxY00 = btYs[qc1]*fracC1+btYs[qc2]*fracC2
    cxZ00 = btZs[qc1]*fracC1+btZs[qc2]*fracC2

    # new centroid of combined cluster
    cxX11 = cX1*fracC1+cX2*fracC2
    cxY11 = cY1*fracC1+cY2*fracC2
    cxZ11 = cZ1*fracC1+cZ2*fracC2

    # offset between two (b/c not all links used)
    dcxX = cxX11-cxX00
    dcxY = cxY11-cxY00
    dcxZ = cxZ11-cxZ00


    # check relative shift of cluster 1 (subtracting possible DC offset)
    # offsets for cluster 1: new loc - old loc
    qX_off1 = cX1 - btXs[qc1]
    qY_off1 = cY1 - btYs[qc1]
    qZ_off1 = cZ1 - btZs[qc1]        
    if abs(qZ_off1-dcxZ) > vshiftmax
        continue
    elseif sqrt((qX_off1-dcxX)^2 + (qY_off1-dcxY)^2) > hshiftmax
        continue
    end 

    # check relative shift of cluster 2 (subtracting possible DC offset)
    # offsets for cluster 2: new loc - old loc
    qX_off2 = cX2 - btXs[qc2]
    qY_off2 = cY2 - btYs[qc2]
    qZ_off2 = cZ2 - btZs[qc2]        
    if abs(qZ_off2-dcxZ) > vshiftmax
        continue
    elseif sqrt((qX_off2-dcxX)^2 + (qY_off2-dcxY)^2) > hshiftmax
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
    brXs[c1qixs] .+= qX_off1
    brYs[c1qixs] .+= qY_off1
    brZs[c1qixs] .+= qZ_off1
    brorgs[c1qixs] .+= qtim_off1

    # update locations: cluster 2
    brXs[c2qixs] .+= qX_off2
    brYs[c2qixs] .+= qY_off2
    brZs[c2qixs] .+= qZ_off2
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
    brXs[iclust] .+= (cxX00 - cxX11)
    btXs[iclust] .= cxX00
    brYs[iclust] .+= (cxY00 - cxY11)
    btYs[iclust] .= cxY00
    brZs[iclust] .+= (cxZ00 - cxZ11)
    btZs[iclust] .= cxZ00
    #btorgs[iclust] .= 0.0 # always zero, never updated
    end

    end # ---- end of clustering loop

    # return results
    bnb = btnbranch[brcids]
    return brXs, brYs, brZs, brorgs, brcids, bnb
end

##################

####### function with i64 ixx arrays
function clustertree_3D(pqix1::Vector{Int32},pqix2::Vector{Int32},ixx1::Vector{Int64},ixx2::Vector{Int64},
    tdif::Vector{Float32},itab::Vector{Int16},
    qXs::Vector{Float64},qYs::Vector{Float64},qZs::Vector{Float64},ttTABs::AbstractVector{},
    nit::Int64,boxwid::Float64,irelonorm::Int64,rmsmax::Float32,rmedmax::Float32,
    distmax::Float64,distmax2::Float64,hshiftmax::Float64,vshiftmax::Float64,torgdifmax::Float32,
    nupdate::Int64,maxlink::Int64)

    # setup parameters
    cZmin = min(0.0,minimum(qZs))
    distmax22 = distmax^2 # in squared km

    # array sizes
    nq = length(qXs)
    bnpair = length(pqix1)

    # initialize relocated event arrays
    brXs = copy(qXs)
    brYs = copy(qYs)
    brZs = copy(qZs)
    brorgs = zeros(Float32,nq)
    brcids = Vector{Int32}(1:nq)

    # initialize clustering tree arrays
    btXs = copy(brXs)
    btYs = copy(brYs)
    btZs = copy(brZs)
    btorgs = zeros(Float32,nq)
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
    cdist22 = (btXs[qc2]-btXs[qc1])^2 + (btYs[qc2]-btYs[qc1])^2 + (btZs[qc2]-btZs[qc1])^2
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
    dqX1, dqY1 = zeros(npick),zeros(npick)
    dqZ1, dqorg1 = zeros(npick),zeros(Float32,npick)
    dqX2, dqY2 = zeros(npick),zeros(npick)
    dqZ2, dqorg2 = zeros(npick),zeros(Float32,npick)
    tab21 = zeros(Int16,npick)
    tdif21 = zeros(Float32,npick)
    ix1 = 1 # start index for event pair in the arrays above
    @inbounds for jj in linx
        pix1, pix2, jx1, jx2 = pqix1[jj],pqix2[jj],ixx1[jj],ixx2[jj]
        ix2 = ix1 + (jx2 - jx1) # end-index in npick array
        if brcids[pix1]==qc1 # regular: event 1 in cluster 1, event 2 in cluster 2
            dqX1[ix1:ix2] .= brXs[pix1]-btXs[qc1]
            dqX2[ix1:ix2] .= brXs[pix2]-btXs[qc2]
            dqY1[ix1:ix2] .= brYs[pix1]-btYs[qc1]
            dqY2[ix1:ix2] .= brYs[pix2]-btYs[qc2]
            dqZ1[ix1:ix2] .= brZs[pix1]-btZs[qc1]
            dqZ2[ix1:ix2] .= brZs[pix2]-btZs[qc2]
            dqorg1[ix1:ix2] .= brorgs[pix1]-btorgs[qc1]
            dqorg2[ix1:ix2] .= brorgs[pix2]-btorgs[qc2]
            tdif21[ix1:ix2] .= tdif[jx1:jx2] # keep the tdif with same sign
        else # flipped: event 1 in cluster 2, event 2 in cluster 1
            dqX1[ix1:ix2] .= brXs[pix2]-btXs[qc1]
            dqX2[ix1:ix2] .= brXs[pix1]-btXs[qc2]
            dqY1[ix1:ix2] .= brYs[pix2]-btYs[qc1]
            dqY2[ix1:ix2] .= brYs[pix1]-btYs[qc2]
            dqZ1[ix1:ix2] .= brZs[pix2]-btZs[qc1]
            dqZ2[ix1:ix2] .= brZs[pix1]-btZs[qc2]
            dqorg1[ix1:ix2] .= brorgs[pix2]-btorgs[qc1]
            dqorg2[ix1:ix2] .= brorgs[pix1]-btorgs[qc2]
            tdif21[ix1:ix2] .= -tdif[jx1:jx2] # flip the tdif
        end
        tab21[ix1:ix2] .= itab[jx1:jx2] # stays same
        ix1 = ix2 + 1 # update start index in npick array
    end

    # unweighted cluster centroid
    cX0 = (btXs[qc1] + btXs[qc2]) / 2.0
    cY0 = (btYs[qc1] + btYs[qc2]) / 2.0
    cZ0 = (btZs[qc1] + btZs[qc2]) / 2.0

    # run difclust (relocation norms 1, 2, 3)
    if irelonorm == 1
        cX1, cY1, cZ1, cX2, cY2, cZ2, 
        cdist, torgdif, resid, rms, rmed, resol = difclust1_3D(
        cX0,cY0,cZ0,tdif21,tab21,
        dqX1, dqY1, dqZ1, dqorg1,
        dqX2, dqY2, dqZ2, dqorg2,
        ttTABs,boxwid,nit)
    elseif irelonorm == 2
        cX1, cY1, cZ1, cX2, cY2, cZ2, 
        cdist, torgdif, resid, rms, rmed, resol = difclust2_3D(
        cX0,cY0,cZ0,tdif21,tab21,
        dqX1, dqY1, dqZ1, dqorg1,
        dqX2, dqY2, dqZ2, dqorg2,
        ttTABs,boxwid,nit)
    else
        cX1, cY1, cZ1, cX2, cY2, cZ2, 
        cdist, torgdif, resid, rms, rmed, resol = difclust3_3D(
        cX0,cY0,cZ0,tdif21,tab21,
        dqX1, dqY1, dqZ1, dqorg1,
        dqX2, dqY2, dqZ2, dqorg2,
        ttTABs,boxwid,nit)
    end

    # careful with cluster mergers near surface
    if (min(cZ1,cZ2) < cZmin)
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
    cxX00 = btXs[qc1]*fracC1+btXs[qc2]*fracC2
    cxY00 = btYs[qc1]*fracC1+btYs[qc2]*fracC2
    cxZ00 = btZs[qc1]*fracC1+btZs[qc2]*fracC2

    # new centroid of combined cluster
    cxX11 = cX1*fracC1+cX2*fracC2
    cxY11 = cY1*fracC1+cY2*fracC2
    cxZ11 = cZ1*fracC1+cZ2*fracC2

    # offset between two (b/c not all links used)
    dcxX = cxX11-cxX00
    dcxY = cxY11-cxY00
    dcxZ = cxZ11-cxZ00


    # check relative shift of cluster 1 (subtracting possible DC offset)
    # offsets for cluster 1: new loc - old loc
    qX_off1 = cX1 - btXs[qc1]
    qY_off1 = cY1 - btYs[qc1]
    qZ_off1 = cZ1 - btZs[qc1]        
    if abs(qZ_off1-dcxZ) > vshiftmax
        continue
    elseif sqrt((qX_off1-dcxX)^2 + (qY_off1-dcxY)^2) > hshiftmax
        continue
    end 

    # check relative shift of cluster 2 (subtracting possible DC offset)
    # offsets for cluster 2: new loc - old loc
    qX_off2 = cX2 - btXs[qc2]
    qY_off2 = cY2 - btYs[qc2]
    qZ_off2 = cZ2 - btZs[qc2]        
    if abs(qZ_off2-dcxZ) > vshiftmax
        continue
    elseif sqrt((qX_off2-dcxX)^2 + (qY_off2-dcxY)^2) > hshiftmax
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
    brXs[c1qixs] .+= qX_off1
    brYs[c1qixs] .+= qY_off1
    brZs[c1qixs] .+= qZ_off1
    brorgs[c1qixs] .+= qtim_off1

    # update locations: cluster 2
    brXs[c2qixs] .+= qX_off2
    brYs[c2qixs] .+= qY_off2
    brZs[c2qixs] .+= qZ_off2
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
    brXs[iclust] .+= (cxX00 - cxX11)
    btXs[iclust] .= cxX00
    brYs[iclust] .+= (cxY00 - cxY11)
    btYs[iclust] .= cxY00
    brZs[iclust] .+= (cxZ00 - cxZ11)
    btZs[iclust] .= cxZ00
    #btorgs[iclust] .= 0.0 # always zero, never updated
    end

    end # ---- end of clustering loop

    # return results
    bnb = btnbranch[brcids]
    return brXs, brYs, brZs, brorgs, brcids, bnb
end