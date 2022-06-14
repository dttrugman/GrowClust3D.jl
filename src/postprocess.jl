########### Function to Finalize Clustering Tree ###########

function make_clustertree(rlons,rlats,rdeps,rorgs,rcids,qix)

    println("\nFinalizing clustering trees...")

    # temporary dataframe with event and cluster number
    tdf = DataFrame("enum"=>qix,"cnum"=>rcids)

    # compute nbranch
    transform!(groupby(tdf, :cnum), nrow => :nb)

    # unique clusters, sort by nbranch
    select!(tdf,[:cnum,:nb])
    unique!(tdf)
    sort!(tdf,[:nb,:cnum],rev=[true,false])

    # assign new cluster ids, starting with largest
    tdf[!,:cid] =range(1,nrow(tdf),step=1)

    # add back in event ids, reorder by event index
    tdf = innerjoin(DataFrame("enum"=>qix,"cnum"=>rcids),tdf,on=:cnum)
    sort!(tdf,:enum)
    nclust = maximum(tdf.cid)

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

    # return results
    return rnbranch, tnbranch,tlats,tlons,tdeps,torgs

end

#################################################################
#################################################################

###### Function to Compute Misfits (w/ otime adjustment) ########

function compute_misfits(inpD,xdf,rdf,ttTABs)
    println("\nComputing misfits...")

    # number of travel time tables
    ntab = length(ttTABs)

    ### Misfits for a 1D Velocity model
    if inpD["tt_ndim"] == 2

        # select columns
        resdf = select(xdf,[:qid1,:qid2,:tdif,:itab,:sX4,:sY4,:rxcor,:sdist])

        # merge with event data, renaming columns to specify 1/2
        resdf = innerjoin(resdf,rdf[!,[:enum,:evid,:rX,:rY,:rdep,:rtim]],on=:qid1=>:evid)
        DataFrames.rename!(resdf,:enum=>:qnum1,:rX=>:qX1,:rY=>:qY1,:rdep=>:qZ1,:rtim=>:qtim1)
        resdf = innerjoin(resdf,rdf[!,[:enum,:evid,:rX,:rY,:rdep,:rtim]],on=:qid2=>:evid)
        DataFrames.rename!(resdf,:enum=>:qnum2,:rX=>:qX2,:rY=>:qY2,:rdep=>:qZ2,:rtim=>:qtim2)

        # here, compute misfits for relocated event pairs and good differential times (for consistency with f90 codes)
        rcids = rdf.rcid
        esdf = resdf[(rcids[resdf.qnum1] .== rcids[resdf.qnum2]) .& ( # only event pairs in same cluster
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

    else # --- misfits for a 3D model

        # select columns
        resdf = select(xdf,[:qid1,:qid2,:tdif,:itab,:rxcor,:sdist])

        # merge with event data, renaming columns to specify 1/2
        resdf = innerjoin(resdf,rdf[!,[:enum,:evid,:rX,:rY,:rdep,:rtim]],on=:qid1=>:evid)
        DataFrames.rename!(resdf,:enum=>:qnum1,:rX=>:qX1,:rY=>:qY1,:rdep=>:qZ1,:rtim=>:qtim1)
        resdf = innerjoin(resdf,rdf[!,[:enum,:evid,:rX,:rY,:rdep,:rtim]],on=:qid2=>:evid)
        DataFrames.rename!(resdf,:enum=>:qnum2,:rX=>:qX2,:rY=>:qY2,:rdep=>:qZ2,:rtim=>:qtim2)

        # here, compute misfits for relocated event pairs and good differential times (for consistency with f90 codes)
        rcids = rdf.rcid
        resdf = resdf[(rcids[resdf.qnum1] .== rcids[resdf.qnum2]) .& ( # only event pairs in same cluster
            resdf[!,:sdist].<=inpD["delmax"]).&(resdf[!,:rxcor].>=inpD["rmin"]),:]

        # compute predicted travel times
        resdf[!,:pdif] .= 0.0
        for ii = 1:nrow(resdf)
            resdf[ii,:pdif] = ttTABs[resdf[ii,:itab]](resdf[ii,:qX2],resdf[ii,:qY2],resdf[ii,:qZ2]) -
                ttTABs[resdf[ii,:itab]](resdf[ii,:qX1],resdf[ii,:qY1],resdf[ii,:qZ1]) +
                resdf[ii,:qtim2] - resdf[ii,:qtim1]
        end

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

    ### now compute event statistics

    # arrays to store stats for each event
    nq = nrow(rdf)
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

    ### return all results
    return npp, nss, rmsP, rmsS, msresP, msresS, 
        qrmsP, qrmsS, qndiffP, qndiffS, qnpair

end

###############################################################
###############################################################

### Compute bootstrap statistics ###

function compute_bootstats(inpD,degkm,rdf,bnbM,blonM,blatM,bdepM,borgM)

    # pre-allocate: defaults are NaN for errors, 0 for nb arrays
    nq = nrow(rdf)
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

    ### return results
    return boot_madH, boot_madZ, boot_madT, boot_stdH, boot_stdZ, boot_stdT,
        boot_nbL, boot_nbM, boot_nbH

end