######## Function to Make 1D Raytrace Tables #########

function make_trace1D_tables(inpD,maxSR,max_selev,usta,sta2elev,ntab,
    erad,vzmodel_type,shallowmode)

    ### In-bounds check
    if maxSR >= inpD["tt_xmax"]
        println("ERROR: station distance exceeds tt_xmax!")
        println("Please modify input file accordingly.")
        exit()
    end

    ### Read in velocity model
    println("\nReading velocity model..")
    # read in
    z_s0, alpha_s0, beta_s0 = read_vzmodel(
        inpD["fin_vzmdl"],vpvs=inpD["vpvs_factor"])
    for ii in 1:length(z_s0) # print out
        @printf("%5.2fkm: %6.4f %6.4f\n",z_s0[ii],alpha_s0[ii],beta_s0[ii])
    end
    # check
    if max_selev > -z_s0[1]
        println("ERROR: station elevation above velocity model start!")
        println("Velocity model start elevation: ", -z_s0[1])
        println("Maximum station elevation:", max_selev)
        exit()
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
    qdeptab = collect(range(inpD["tt_zmin"],inpD["tt_zmax"],step=inpD["tt_zstep"]))
    sdeltab = collect(range(inpD["tt_xmin"],inpD["tt_xmax"],step=inpD["tt_xstep"]))

    # store all tables here
    ttLIST = [] # stations for P, then stations for S
    phases = ["P","S"] # make P and S tables
    plongcuts = [inpD["plongcutP"],inpD["plongcutS"]] # cutoffs

    # loop over phases
    println("STARTING RAY TRACING:")
    total_time = @elapsed for (iphase, phase) in enumerate(phases)

        # Print results
        println("\nWorking on Phase: $phase")

        # loop over stations 
        for ista = 1:Int64(ntab/2)

            # get station depth (negative of elevation)
            sta = usta[ista] # station name
            zstart = -sta2elev[sta] # elevation
            println("Working on phase $phase station $sta: depth $zstart km")
            ttabfile = @sprintf("tt.%s.%sg",sta,lowercase(phase))
            
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
                write_smtrace_table(inpD["fdir_ttab"] * ttabfile,inpD["fin_vzmdl"],iphase,
                vzmodel_type,TT,qdeptab, sdeltab,ptab,zstart)
            end

        end # end loop over stations
        
    end
    
    # done
    @printf("\nElapsed seconds: %.2f",total_time)
    println()

    # return final list
    return ttLIST

end

#############################################################
#############################################################

######## Function to Tables from 1D or 3D NLL Grids #########

function make_nllgrid_tables(inpD,maxSR,usta,shallowmode)

    # store all tables here
    nstaU = length(usta)
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

            # Alternate version (if needed)
            if phase == "S"
                aphase = "P"
                vscale = inpD["vpvs_factor"] # S grid = P grid x vscale
            else
                aphase = "S"
                vscale = 1.0/inpD["vpvs_factor"] # P grid = S grid x vscale
            end
            aname = inpD["fin_vzmdl"] * "." * aphase * "." * sta *".time"

            # get header
            if isfile(inpD["fdir_ttab"]* fname * ".hdr")
                grdparams = read_nll_head(inpD["fdir_ttab"]*fname)
                grdparams["vscale"] = Float32(1.0)
            elseif isfile(inpD["fdir_ttab"]* aname * ".hdr")
                fname = aname # update file name to alternate
                grdparams = read_nll_head(inpD["fdir_ttab"]*fname)
                grdparams["vscale"] = Float32(vscale) # apply Vp/Vs to existing grid
            else
                println("ERROR: MISSING GRID: $fname")
                exit()
            end

            # check projection
            proj_ok = check_proj(grdparams, inpD)
            if !(proj_ok)
                println("ERROR: PROJECTION MISMATCH: $fname")
                exit()
            end

            # check bounds (only needed for TIME2D)
            if grdparams["gtype"] == "TIME2D"
                gmaxR = grdparams["yORG"] + grdparams["dY"]*Float64(grdparams["nY"]-1)
                if gmaxR < maxSR
                    println("ERROR: maximum station distance too large for grid!")
                    println("max station distance: $maxSR, max grid distance: $gmaxR")
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

    # done
    @printf("\nElapsed seconds: %.2f",total_time)
    println()

    # return final list
    return ttLIST
end

