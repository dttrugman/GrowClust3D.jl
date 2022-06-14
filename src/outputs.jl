
######### Function to Write Catalog File ###########

function write_cat(inpD,rdf,qnpair,qndiffP,qndiffS,
    qrmsP,qrmsS,boot_madH,boot_madZ,boot_madT)

    println("\nWriting output catalog: ", inpD["fout_cat"])

    # open output file
    fcat = open(inpD["fout_cat"],"w")

    # loop over events
    for ii = 1:nrow(rdf)
        
        # print out origin time, relocated position, magnitude
        dateS = Dates.format(rdf[ii,:rot],"YYYY mm dd HH MM SS.sss")
        @printf(fcat,"%s %9d %9.5f %10.5f %7.3f %5.2f ",
            dateS,rdf[ii,:evid],rdf[ii,:rlat],rdf[ii,:rlon],
            rdf[ii,:rdep],rdf[ii,:mag])
        
        # print out cluster number and fits
        @printf(fcat,"%7d %7d %7d %5d %5d %5d %5.2f %5.2f ",
            rdf[ii,:enum],rdf[ii,:rcid],rdf[ii,:rnb],qnpair[ii],
            qndiffP[ii],qndiffS[ii],qrmsP[ii],qrmsS[ii])
        
        # print out uncertanties and catalog locations
        @printf(fcat,"%7.3f %7.3f %7.3f %9.5f %10.5f %7.3f\n",
            boot_madH[ii],boot_madZ[ii],boot_madT[ii],
            rdf[ii,:qlat],rdf[ii,:qlon],rdf[ii,:qdep])
    end

    # close file
    close(fcat)

end

####################################################
####################################################

######### Function to Write Cluster File ###########

function write_clust(inpD,rdf,tnbranch,tlats,tlons,tdeps,torgs,
    boot_madH, boot_madT, boot_madZ)

    println("\nWriting output clusters")

    # open output file
    fcc = open(inpD["fout_clust"],"w")

    # loop over all clusters to output (only makes sense for n>=2)
    for cc in findall(tnbranch .>= min(2,inpD["nbranch_min"]))
        
        # write cluster info
        @printf(fcc,"%8d %7d %9.5f %10.5f %7.3f %7.3f\n",
            cc,tnbranch[cc],tlats[cc],tlons[cc],tdeps[cc],torgs[cc])
        
        # write info for all events
        for ii in findall(rdf.rcid==cc)
        
            # print out clustering, event info, mag, otime
            dateS = Dates.format(rdf[ii,:rot],"YYYY mm dd HH MM SS.sss")
            @printf(fcc,"%8d %8d %9d %5.2f %s ",cc,ii,rdf[ii,:evid],rdf[ii,:mag],dateS)
            
            # print event location and position w/in cluster
            qdx = degkm*(rdf[ii,:rlon]-tlons[cc])*cosd(tlats[cc])
            qdy = degkm*(rdf[ii,:rlat]-tlats[cc])
            qdz = rdf[ii,:rdep]-tdeps[cc]
            @printf(fcc,"%9.5f %10.5f %7.3f %9.4f %9.4f %9.4f ",
                rdf[ii,:rlat],rdf[ii,:rlon],rdf[ii,:rdep],qdx,qdy,qdz)

            # print out uncertanties  and catalog locations
            @printf(fcc,"%7.3f %7.3f %7.3f %9.5f %10.5f %7.3f\n",
                boot_madH[ii],boot_madZ[ii],boot_madT[ii],
                rdf[ii,:qlat],rdf[ii,:qlon],rdf[ii,:qdep])
        end
        
    end

    # close file 
    close(fcc)

end

####################################################
####################################################

######### Function to Write Bootstrap File #########

function write_boot(inpD,rdf,boot_madH,boot_madT,boot_madZ,
    boot_nbH,boot_nbL,boot_nbM,boot_stdH,boot_stdT,boot_stdZ,
    blatM,blonM,bdepM,bnbM)
    
    println("\nWriting output bootstrapping")
    
    # open file
    fbb = open(inpD["fout_boot"],"w")
    
    # file header
    nq = nrow(rdf)
    @printf(fbb,"%8d %5d\n",nq,inpD["nboot"])
    
    # loop over events
    for ii = 1:nq
        
        # event header
        @printf(fbb,"%8d %9d %9.5f %10.5f %7.3f ",
            rdf[ii,:enum],rdf[ii,:evid],rdf[ii,:rlat],rdf[ii,:rlon],rdf[ii,:rdep])
        @printf(fbb,"%7d %7d %6.2f %7d %7d %9.5f %10.5f %7.3f ",
            rdf[ii,:rcid],rdf[ii,:rnb],boot_nbM[ii],boot_nbL[ii],boot_nbH[ii],
            rdf[ii,:qlat],rdf[ii,:qlon],rdf[ii,:qdep])
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

####################################################
####################################################

########## Function to Write Log File ##############

function write_log(inpD,auxparams,runstats,tnbranch)

    println("\nWriting run log")

    # unpack aux parameters
    infile_ctl, distmax, distmax2, hshiftmax, vshiftmax, rmedmax = auxparams

    # unpack run statistics
    nq, nreloc, npair, qnpair, npp, nss, rmsP, rmsS, msresP, msresS = runstats

    # compute tree counts
    ntree2 = sum(tnbranch.>=2)
    ntree5 = sum(tnbranch.>=5)
    ntree10 = sum(tnbranch.>=10)
    ntree20 = sum(tnbranch.>=20)
    ntree50 = sum(tnbranch.>=50)
    ntree100 = sum(tnbranch.>=100)

    # write log to output file
    if !(inpD["fout_log"] in ["none","None", "NONE"])
        
        # open file
        flog = open(inpD["fout_log"],"w")

        # Run Parameters
        @printf(flog, "********************** Input Parameters ***********************\n")
        @printf(flog, "     control file:   %s\n", infile_ctl)
        @printf(flog, "       event list:   %s\n", inpD["fin_evlist"])
        @printf(flog, "     station list:   %s\n", inpD["fin_stlist"])
        @printf(flog, "     xcordat file:   %s\n", inpD["fin_xcordat"])
        @printf(flog, "  travel time src:   %s\n", inpD["ttabsrc"])
        @printf(flog, "     velocity mdl:   %s\n", inpD["fin_vzmdl"])
        @printf(flog, "   map projection:   %s %s %.6f %.6f %.2f\n", 
            inpD["proj"], inpD["rellipse"], inpD["lon0"], inpD["lat0"], inpD["rotANG"])
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

    # otherwise, write report to terminal
    else

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
        @printf("************************ Output files *************************\n")
        @printf("     catalog file:   %s\n", inpD["fout_cat"])
        @printf("     cluster file:   %s\n", inpD["fout_clust"])
        @printf("         log file:   %s\n", inpD["fout_log"])
        @printf("   bootstrap file:   %s\n", inpD["fout_boot"])
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
        @printf("\n")
        @printf("==================================================================\n")
        @printf("==================================================================\n")
        @printf("\n")
        @printf("********************  GROWCLUST Run Summary  *********************\n")
        @printf("%55s %10d\n", "Number of catalog events: ", nq)
        @printf("%55s %10d\n", "Number of relocated events: ", nreloc)
        @printf("%55s %10d\n", "Number of input event pairs: ", npair)
        @printf("%55s %10d\n", "Number of event pairs used: ", sum(qnpair)/2)
        @printf("%55s %10d\n", "Number of xcor data used (total, P+S): ", npp + nss)
        @printf("%55s %10d\n", "Number of xcor data used (P-phase): ", npp)
        @printf("%55s %10.4f\n", "RMS differential time residual (P-phase): ", rmsP)
        @printf("%55s %10.4f\n", "Mean (signed) differential time residual (P-phase): ", msresP)
        @printf("%55s %10d\n", "Number of xcor data used (S-phase): ", nss)
        @printf("%55s %10.4f\n", "RMS differential time residual (S-phase): ", rmsS)
        @printf("%55s %10.4f\n", "Mean (signed) differential time residual (S-phase): ", msresP)
        @printf("\n")
        @printf("%55s %9d\n", "Number of clusters with >=   2 events: ", ntree2)
        @printf("%55s %9d\n", "Number of clusters with >=   5 events: ", ntree5)
        @printf("%55s %9d\n", "Number of clusters with >=  10 events: ", ntree10)
        @printf("%55s %9d\n", "Number of clusters with >=  20 events: ", ntree20)
        @printf("%55s %9d\n", "Number of clusters with >=  50 events: ", ntree50)
        @printf("%55s %9d\n", "Number of clusters with >= 100 events: ", ntree100)

    end
end