### READ_VZMODEL reads velocity model from text file
#   The format of file is Z, Vp, Vs at each depth point.
#   Vs is optional, and can be instead computed from Vp and vpvs_factor
#   The assembled velocity model is meant to be interpolated, so 
#   layer interfaces need to be explicitly defined with two depth points
#   with top and bottom layer velocities.
#
#  > Inputs
#     fname:       Name of velocity model file
#     vpvs: Assumed Vp/Vs to use if Vs is not listed
#
#  < Returns
#     zz: Depth points
#     vp: Vp at each depth
#     vs: Vs at each depth
#     const_VpVs: true or false, if Vs is just a constant times Vp

function read_vzmodel(fname; vpvs=sqrt(3.0))
   
    # read all lines in file
    lines = readlines(fname)
    npts = length(lines)
    
    # initialize arrays
    zz = Array{Float64}(undef,npts)
    alpha = Array{Float64}(undef,npts)
    beta = Array{Float64}(undef,npts)
    
    # parse file
    for (ii,line) in enumerate(lines)
       z, p, s = split(line)
       zz[ii] = parse(Float64,z)
       alpha[ii] = parse(Float64,p)
       b = parse(Float64,s)
       if b == 0.0
           b = alpha[ii]/vpvs
       end
       beta[ii] = b
    end

    # check to see if constant vp/vs
    beta0 = alpha ./ vpvs
    if maximum(abs.(beta-beta0)) < 0.001
        const_VpVs = true
    else
        const_VpVs = false
    end
    
    # return data
    return zz, alpha, beta, const_VpVs
    
end

#######################################################################


### FIND_MOHO: Function to find moho within velocity model
#
#  > Inputs: (zz, alpha), depth points and Vp in velocity model
#  < Returns: imoho, index in velocity model of Moho depth (0 = not found)

function find_moho(zz, alpha; alphaM1=7.2, alphaM2=8.0, dalphaM=0.5)
    
    # initialize
    imoho = findfirst(alpha .>= 0.5*(alphaM1 + alphaM2)) # set a reasonable default
    npts = length(zz)

    # now look for jump in Vp, above ~ 7.2km/s and at most 8.0 km/s (by default)
    #  (this happens at >~ 30km in continental crust, but much shallower in oceans)
    for ii in 2:npts
        if ((alpha[ii] >= alphaM1) & (alpha[ii-1] < alpha[ii]-dalphaM))
            imoho=ii
            break
        elseif (alpha[ii] >= alphaM2)
            imoho=ii
            break
        end
    end
    
   # return index
   return imoho
end

### INTERP_VZMODEL creates a finely spaced version while preserving layer interfaces
#   The input velocity model is meant to be interpolated, so 
#   layer interfaces need are explicitly defined with two depth points
#   with top and bottom layer velocities.
#
#  > Inputs
#     z_s0:      Depth points in coarse model
#     alpha_s0:  P-wave speed in coarse model
#     beta_s0:   S-wave speed in coarse model
#     itp_dz:    Depth interval for interpolation
#     ztop:      Top depth of the output model (can be higher)
#
#  < Returns
#     z_s:       Depth points in finer interpolated model
#     alpha_s:   P-wave speed in finer interpolated model
#     beta_s:    S-wave speed in finer interpolated model

function interp_vzmodel(z_s0, alpha_s0, beta_s0; itp_dz=1.0, ztop=-3.0)
    
    # define evenly spaced grid spanning velocity model
    zmin = min(floor(ztop),floor(z_s0[1]),-3.0) # down to -3.0 for a pad
    zmax = ceil(z_s0[end])
    zitp = range(zmin, zmax, step=itp_dz)

    # output z points: points from grid while preserving known interfaces
    z_s = sort(vcat(setdiff(zitp,z_s0), z_s0))
    ns = length(z_s)
    alpha_s, beta_s = zeros(Float64,ns), zeros(Float64,ns)

    # interpolate to output points while preserving interfaces
    ii=1 
    ns0 = length(z_s0)
    for jj = 1:length(z_s)
        if z_s[jj] < z_s0[1] # top extrapolate
            alpha_s[jj]=alpha_s0[1]
            beta_s[jj]=beta_s0[1]
        elseif (ii>ns0) # bottom extrapolate
            alpha_s[jj]=alpha_s0[end]
            beta_s[jj]=beta_s0[end]
        elseif z_s[jj]==z_s0[ii] # preserve interface
            alpha_s[jj]=alpha_s0[ii]
            beta_s[jj]=beta_s0[ii]
            ii+=1
        else # linear interpolation for new grid
            zfrac = (z_s[jj] - z_s0[ii-1])/(z_s0[ii]-z_s0[ii-1])
            alpha_s[jj] = alpha_s0[ii-1] + zfrac*(alpha_s0[ii]-alpha_s0[ii-1])
            beta_s[jj] = beta_s0[ii-1] + zfrac*(beta_s0[ii]-beta_s0[ii-1])
        end
    end

    # # dummy interface at the bottom???
    # push!(z_s, z_s[end])
    # push!(alpha_s, alpha_s[end])
    # push!(beta_s, beta_s[end])
    
    # return
    return z_s, alpha_s,beta_s 
end

#######################################################################


### EFLATTEN calculates flat earth tranformation
# > Inputs: (z_s,vel_s) in stardard (spherical earth) coords.
# < Returns: (z_f,vel_f) in flat-earth coords.
function eflatten(z_s, vel_s; erad = 6371.0)
    r = erad .- z_s
    z_f = -erad .* log.(r/erad)
    vel_f = vel_s .* (erad./r)
    return z_f, vel_f
end

### UNFLATTEN is inverse of FLATTEN.
#  > Inputs: (z_f,vel_f) in flat-earth coords.
#  < Return: (z_s,vel_s) in stardard (spherical earth) coords.
function uneflatten(z_f, vel_f; erad=6371.0)
    r = erad*exp.(-z_f./erad)
    z_s = erad.-r
    vel_s = vel_f.*(r./erad)
    return z_s, vel_s
end

#######################################################################


### LAYER_TRACE calculates the travel time and range offset 
###   for ray tracing through a single layer.
#
#  > Inputs  
#     p:      horizontal slowness
#     h:     layer thickness
#     utop:  slowness at top of layer
#     ubot:  slowness at bottom of layer
#     imth:  interpolation method for integration
#              imth = 1,  v(z) = 1/sqrt(a - 2*b*z)     fastest to compute
#                   = 2,  v(z) = a - b*z               linear gradient
#                   = 3,  v(z) = a*exp(-b*z)           preferred when Earth Flattening is applied
#  < Returns  
#     dx:   range offset
#     dt:   travel time
#     irtr: return code
#           = -1, zero thickness layer
#           =  0,  ray turned above layer
#           =  1,  ray passed through layer
#           =  2,  ray turned within layer, 1 segment counted
#
#
### Ray tracing background ###
#
#  Total slowness u = 1 / v
#  Horizontal slowness = p = u(z)*sin(theta) = dT/dX = u_tp = constant for a given ray
#  Vertical slowness = eta = sqrt(u**2-p**2)
#  Theta = ray incidence from vertical (0 = straight down, pi/2=90deg = horizontal)
# 
#  For each layer, integrate ray-tracing equations:
#   dT = integral(u^2/eta,dz)
#   dX = integral(p/eta,dz)
#   dTau = integral(eta,dz)
#    where Tau(p) = T(p) - p*X(p)
#
function layer_trace(p, h, utop, ubot, imth)
    
    # check for zero thickness layer
    if (h == 0.0)
        dx=0.0
        dt=0.0
        irtr=-1
        return dx, dt, irtr
    end
    
    # slowness at top of layer
    u = utop
    
    # check for complex vertical slowness: ray turned above layer
    if (u-p <= 0.0)
        dx=0.0
        dt=0.0
        irtr = 0 # return code: ray turned above layer
        return dx, dt, irtr
    end
    
    # calculate vertical slowness, eta = sqrt(u^2-p^2)
    eta2=(u-p)*(u+p)
    eta = sqrt(eta2)
    
    # special function needed for integral at top of layer
    if (imth == 2)
        y=u+eta
        if p != 0.0; y=y/p; end
        qr=log(y)
    elseif (imth == 3)
        qr=atan(eta,p) # inverse tangent of eta/p
    end 
        
    # b factor (ray tracing integral constant in denominator)
    if (imth == 1)
        b=-(utop^2-ubot^2)/(2.0*h)
    elseif (imth == 2)         
        vtop=1.0/utop
        vbot=1.0/ubot
        b=-(vtop-vbot)/h
    else                     
        b=-log(ubot/utop)/h   # flat earth     
    end
    
    # constant velocity layer: no need to integrate - just multiply
    if (b == 0.0)
        b=1.0/h          # b = 1/dz
        etau=eta        # integrand for Tau equation
        ex=p/eta        # integrand for X equation
        irtr=1          # return code: ray passed through layer
    
    # non-constant layer: need to integrate
    else    
    
        # ray tracing integral at upper limit, 1/b factor omitted until end
        if (imth == 1)
            etau=-eta2*eta/3.0
            ex=-eta*p
        elseif (imth == 2)
            ex=eta/u
            etau=qr-ex
            if p != 0.0; ex=ex/p; end
        else
            etau=eta-p*qr
            ex=qr
        end

        # check lower limit for turning point
        u = ubot
        if (u <= p)
            irtr=2          # return code: ray turned within layer
            
        # no turning point: ray passes through 
        else
            irtr=1            # return code: ray passed through layer
            eta2=(u-p)*(u+p)  # update vertical slowness for ubot
            eta=sqrt(eta2) # update vertical slowness for ubot
            if (imth == 1)
                etau=etau+eta2*eta/3.0 # Tau equation
                ex=ex+eta*p           # X equation
            elseif (imth == 2)
                y=u+eta
                z=eta/u
                etau=etau+z
                if (p != 0.0)
                    y=y/p
                    z=z/p
                end
                qr=log(y)
                etau=etau-qr        # Tau equation
                ex=ex-z             # X equation
            else
                qr=atan(eta,p)
                etau=etau-eta+p*qr  # Tau equation
                ex=ex-qr            # X equation
            end
        end
    end

    # finish ray-tracing equations to get dx, dt
    dx=ex/b         # horizontal offset in km
    dtau=etau/b     # delay time
    dt=dtau+p*dx    # convert delay time to travel time
    return dx, dt, irtr # return
end

#######################################################################


### TRACE_RAYS: Trace rays with different ray parameter p, 
#    tracking offsets and travel times to different depths.
#
#  > Inputs
#     iw:       Phase index in slowness array (0=P, 1=S)
#     z_s:      Depth points in interpolated velocity model
#     z:        Depth points in flattened model
#     slow:     Array of slownesses at each depth point
#     qdeptab:  Array of source depths in output table
#     itp_dz:   Interpolation depth spacing
#     zstart:   Starting depth for "surface" (i.e negative of selev)
#
#  < Returns
#     ptab:     Array of ray parameters used in ray tracing
#     qdepxcor: Matrix (nray x ndep)of offsets to each source depth, for each ray 
#     qdeptcor: Matrix (nray x ndep) of travel times to each soruce depth, for each ray
#     qdepucor: Matrix (nray x ndep) of slownesses at each source depth, for each ray
#     del2W:   Surface-to-surface offset for each ray
#     tt2W:    Surface-to-surface travel time for each ray
#
function trace_rays(iw,z_s,z,slow,qdeptab,itp_dz,zstart)
    
    # validate zstart
    if zstart < z_s[1]
        println("Error in trace_rays: station above vzmodel start")
        println("$zstart < ",z_s[1])
        exit()
    end

    # preferred interpolation method for earth-flattening
    imth = 3
    erad = 6371.0

    # array sizes
    npts = length(z_s)                 # number of depth points in velocity models
    zmax = 9999.0                      # maximum depth point, used only for robustness
    ndep = length(qdeptab)             # source depth points in travel time table
    nray = 10001                       # number of rays to shoot
    
    # setup arrays
    pmin, pmax = 0.0, slow[1,iw]       # min and max ray parameter
    ptab = collect(range(
            pmin,pmax,length=nray))    # rays to shoot
    qdepxcor = fill(NaN64,(nray,ndep)) # offset x for each (p,z)
    qdeptcor = fill(NaN64,(nray,ndep)) # time t for each (p,z)
    qdepucor = fill(NaN64,(nray,ndep)) # slowness u for each (p,z)
    del2W = fill(NaN64,nray)           # surface-to-surface offset X for each p
    tt2W = fill(NaN64,nray)            # surface-to-surface travel time T for each p
    
    # surface correction (no longer needed?)
    ifix = qdeptab.==zstart
    qdepxcor[:,ifix].=0.0
    qdeptcor[:,ifix].=0.0
    qdepucor[:,ifix].=slow[1,iw]

    # loop over ray parameters
    for (ip, p) in enumerate(ptab)

        # start at 0,0
        x, t = 0., 0.

        # find final depth point above the station
        istart = findlast(y->y<=zstart,z_s)

        # loop over layers
        for ii = istart:(npts-1)

            # check for z>zmax (for robustness)
            if z_s[ii] >= zmax
                println("Error: z_s[ii] > zmax")
                exit()
            end

            # layer thickness
            if ii == istart
                h = z[ii+1] - zstart # may be partway through layer
                if (z[ii+1]==z[ii])
                    slowii = slow[ii,iw]
                else
                    slowii = slow[ii,iw] + h/(z[ii+1]-z[ii])*(slow[ii+1,iw]-slow[ii,iw])
                end 
            else
                h = z[ii+1] - z[ii] # full layer thickness
                slowii = slow[ii,iw]
            end
            if h <= 0.0; continue; end # skip if interface

            # ray tracing through this layer
            dx, dt, irtr = layer_trace(p,h,slowii,slow[ii+1,iw],imth)

            # update x, t after tracing through layer
            x+=dx
            t+=dt

            # exit when ray has turned
            if ((irtr==0) | (irtr==2)); break; end

            # save current x,t,u for ray sampling source depths (stored in deptab)
            idep = findall(abs.(qdeptab.-z_s[ii+1]) .<= itp_dz/1.99) # update for robustness
            if sum(idep)>0
                qdepxcor[ip,idep].=x
                qdeptcor[ip,idep].=t
                qdepucor[ip,idep].=slow[ii+1,iw] # since the ray hasn't turned, we are in layer ii+1
            end
        end
        #----- end loop on layers------#

        # compute final surface-to-surface two-way offset and travel times for this ray
        del2W[ip]=2.0*x
        tt2W[ip]=2.0*t
    end

    # return results
    return ptab, qdepxcor, qdeptcor, qdepucor, del2W, tt2W
    
end

#######################################################################


#### FIRST_ARRIVALS: Function to make table of first arrivals, given ray-tracing resutls
#
#  > Inputs
#     itype:    Travel time table type: 1: (x=km,t=sec), 2=(x=deg,t=min)
#     plongcut: Minimum ray parameter, optionally used to filter out refracted rays
#     qdeptab:  Array of source depths in output table
#     sdeltab:  Array of station distances in output table
#     usurf:    Slowness at the surface of the velocity model
#     zstart:   Station depth for start of ray tracing (i.e negative of selev)
#     ptab:     Array of ray parameters used in ray tracing
#     qdepxcor: Matrix (nray x ndep) of offsets to each source depth, for each ray
#     qdeptcor: Matrix (nray x ndep) of travel times to each soruce depth, for each ray
#     qdepucor: Matrix (nray x ndep) of slownesses at each source depth, for each ray
#     del2W:    Surface-to-surface offset for each ray
#     tt2W:     Surface-to-surface travel time for each ray
#
#  < Returns
#     tt:       Travel time table (ndel,ndep)
#     aa:       Takeoff angle table (ndel, dep) in degrees from vertical (0 = straight down)
#
function first_arrivals(itype, plongcut, qdeptab, sdeltab, usurf, zstart,
                   ptab, qdepxcor, qdeptcor, qdepucor, del2W, tt2W)
    
    # earth / geographic parameters
    erad = 6371.0
    ecircum = 2.0*pi*erad
    kmdeg=ecircum/360.0
    
    # array sizes
    ndep, ndel, nray = length(qdeptab), length(sdeltab), length(ptab)
    
    # track x,t,p,u values for each depth
    ncount = 2*nray # accounts for upgoing and downgoing for each ray parameter
    xsave, tsave = zeros(ncount), zeros(ncount)
    psave, usave = zeros(ncount), zeros(ncount)

    # initialize travel time and takeoff angle tables (X x Z)
    tt = fill(NaN64,(ndel,ndep))
    aa = fill(NaN64,(ndel,ndep))
    
    # loop over depth ranges
    for idep in 1:ndep
        #@printf("Table depth %d/%d\n",idep,ndep)

        icount = 1 # reset index in save array
        xold = 0.0 # current x value

        # upgoing rays from source
        if (qdeptab[idep] > zstart) # skip for sources at zero depth
            for ip in 1:nray
                if isnan(qdepxcor[ip,idep]); break; end
                if isnan(del2W[ip]); break; end 
                if (qdepxcor[ip,idep]<xold); break; end # stop when ray turns inward
                xsave[icount]= qdepxcor[ip,idep] # x = surface-to-source
                tsave[icount]= qdeptcor[ip,idep] # t = surface-to-source
                usave[icount]= qdepucor[ip,idep] # u = slowness at source depth
                psave[icount]= -ptab[ip] # upgoing rays from source are flagged as negative here
                xold = xsave[icount]
                icount+=1 # increment counter
            end
            ip2 = icount # start next loop here
        else
            ip2 = nray # start next loop here
        end

        # downgoing rays from source
        for ip in range(ip2,1,step=-1)
            if isnan(qdepxcor[ip,idep]); continue; end
            if isnan(del2W[ip]); continue; end
            xsave[icount]= del2W[ip] - qdepxcor[ip,idep] # x = (surface-to-surface) minus (surface-to-source)
            tsave[icount]= tt2W[ip] - qdeptcor[ip,idep] #  t = (surface-to-surface) minus (surface-to-source)
            usave[icount]= qdepucor[ip,idep] # u = slowness at source depth
            psave[icount]= ptab[ip] # downgoing rays from source are flagged as positive here
            icount+=1 # increment save counter
        end
        ncount = icount-1
        #println(ncount)

        # loop over offsets
        for (idel, sdel0) in enumerate(sdeltab)
            
            # convert from degrees to km (xsave is always in km)
            if itype == 2 
                sdel = sdel0*kmdeg
            else
                sdel = sdel0
            end

            # search for first-arriving ray at this offset
            tt[idel,idep] = 9999.0
            for ii in 2:ncount

                # get x-span, skipping unless we are bracketing desired offset
                x1, x2 = xsave[ii-1], xsave[ii] # these are always in km
                if ((x1>sdel) | (x2<sdel)); continue; end

                # skip downgoing rays (p > 0) with this incidence: from refracted rays
                if ((psave[ii] > 0.0) & (psave[ii] < plongcut)); continue; end

                # compute tbest
                frac=(sdel-x1)/(x2-x1)
                tbest=tsave[ii-1]+frac*(tsave[ii]-tsave[ii-1])
                
                # check to see if this is the first arrival so far
                if tbest < tt[idel,idep]
                    
                    # update travel time table
                    tt[idel,idep] = tbest
                    
                    # compute pbest and ubest, careful around p = 0
                    if (psave[ii-1]*psave[ii] > 0.0) # both positive or both negative
                    #if (((psave[ii-1]>0.0) & (psave[ii]>0.0)) | ((psave[ii-1]<0.0) & (psave[ii]<0.0)))
                        pbest = psave[ii-1]+frac*(psave[ii]-psave[ii-1])
                        ubest = usave[ii-1]+frac*(usave[ii]-usave[ii-1])
                    elseif (frac < 0.5)
                        pbest = psave[ii-1]
                        ubest = usave[ii-1]
                    else
                        pbest = psave[ii]
                        ubest = usave[ii]
                    end
                        
                    # compute takeoff angle at the source (straight down = 0)
                    angle = asind(abs(pbest)/ubest) # arcsin in degrees
                    if pbest < 0.0; angle = 180.0-angle; end # upgoing ray, need to flip
                    aa[idel,idep] = angle
                
                end
            end
                    
            # no ray arrivals
            if tt[idel,idep] == 9999.0
                tt[idel,idep] = NaN64
                aa[idel,idep] = NaN64
            end
        end
    end
                
    # fix cases for qdep < station
    for idep = 1:ndep
        if qdeptab[idep] > zstart 
            break # stop when we depth > station
        else
            dzdep = zstart - qdeptab[idep] # depth difference
            if itype == 1 # flat earth case
                tt_surf = usurf*sqrt.(dzdep^2 .+ sdeltab.^2) # straight ray
                aa_surf = 90.0 .+ atand.(dzdep,sdeltab) # straight ray angle
            else # spherical earth case
                xx_surf = sdeltab*kmdeg # offset in km
                tt_surf = usurf*sqrt.(dzdep^2 .+ xx_surf.^2) # straight ray tt
                aa_surf = 90.0 .+ atand.(dzdep,xx_surf) # straight ray angle
            end
            ifix = ((isnan.(tt[:,idep])) .| (tt[:,idep] .>= tt_surf))
            tt[:,idep] .= ifelse.(ifix,tt_surf,tt[:,idep]) # update tt
            aa[:,idep] .= ifelse.(ifix,aa_surf,aa[:,idep]) # update aa
        end
    end

        
    # fix edge case: straight up at zero range
    if abs(sdeltab[1])<0.01
        aa[1,:] .= 180.0
    end
            
    # spherical earth
    if itype == 2
        tt ./= 60.0 # to minutes
        aa .*= erad # spherical ray parameter
    end
        
    # return results
    return tt, aa
    
end

#######################################################################


### WRITE_SMTRACE_TABLE: Function to output travel time table to a file
#
#  > Inputs
#     outfile:  Full path to output file
#     vmodel:   Name of velocity model
#     iw:       Phase index in slowness array (0=P, 1=S)
#     itype:    Travel time table type: 1: (x=km,t=sec), 2=(x=deg,t=min)
#     TT:       Travel time table (ndel, ndep)
#     qdeptab:  Array of source depths in output table
#     sdeltab:  Array of station distances in output table
#     ptab:     Array of ray parameters used in ray tracing
#     zstart:   Station depth where rays started (i.e negative of selev)
#
#  < Returns: NONE

function write_smtrace_table(outfile,vmodel,iw,itype,TT,qdeptab,sdeltab,ptab,zstart)
    
    # array sizes
    ndel, ndep = size(TT)
    @assert length(qdeptab) == ndep
    @assert length(sdeltab) == ndel
    
    # open file
    f = open(outfile, "w")
    
    # write header
    @printf(f,"From deptable, file= %20s iw =%2d pmin=%8.5f pmax=%8.5f np=%6d selev=%.3f\n",
        vmodel,iw,ptab[1],ptab[end],length(ptab),-zstart)

    # write table size
    @printf(f,"%5d%5d\n",ndel,ndep)

    # write source depths
    line = @printf(f,"%9s","")
    for dep in qdeptab
        @printf(f,"%9.1f",dep)
    end
    @printf(f,"\n")

    # write travel times at each offset
    for (i, sdel) in enumerate(sdeltab)

        # write offset
        if itype == 1
            @printf(f,"%9.4f",sdel) # in km
        else
            @printf(f,"%9.3f",sdel) # in deg
        end    
        
        # write travel times
        for j in 1:ndep
            @printf(f,"%9.4f",TT[i,j])
        end
        
        # end line
        @printf(f,"\n")
        
    end
                    
    # close file
    close(f)
    
    
end

#######################################################################

### read_smtrace_table: Function to read SeismoTrace Table from a text file

#  > Inputs
#     fname:    Full path to travel time table
#     shallowmode: Sets the extrapolation condition at shallow depth
#     dtype: Output data type for travel time table
#
#  < Returns: NONE
#     TTinterp: interpolation objection, used like tt = TTinterp(sdist,qdep)

function read_smtrace_table(fname,shallowmode="throw",dtype=Float64)
    
    # read all lines in file
    lines = readlines(fname)
    
    # travel time table setup (line 2, skips header)
    sline = split(lines[2])
    ndel = parse(Int64,sline[1])
    ndep = parse(Int64,sline[2])
    dists = zeros(ndel)
    tt = zeros(dtype,(ndel,ndep))
    idel = 1
    
    # source depth setup
    sline = split(lines[3])
    depths = collect([parse(Float64,xx) for xx in sline])
    
    # loop over travel time lines
    for idel in 1:ndel
        sline = split(lines[3+idel])
        dists[idel] = parse(Float64,sline[1])
        tt[idel,:] = collect([parse(dtype,xx) for xx in sline[2:end]])
    end
    
    # define interpolant object, with extrapolation options
    if shallowmode == "reflect"
        TTinterp = LinearInterpolation((dists,depths),tt,
            extrapolation_bc=((Line(),Line()),(Reflect(),Line())))
    elseif shallowmode == "flat"
        TTinterp = LinearInterpolation((dists,depths),tt,
            extrapolation_bc=((Line(),Line()),(Flat(),Line())))
    elseif shallowmode == "line"
        TTinterp = LinearInterpolation((dists,depths),tt,
            extrapolation_bc=((Line(),Line()),(Line(),Line())))
    else
        TTinterp = LinearInterpolation((dists,depths),tt,
            extrapolation_bc=((Line(),Line()),(Throw(),Line())))
    end
    
end

### make_smtrace_table: Function to make SeismoTrace Table from an Array

#  > Inputs
#     dists, depths: 1D arrays of distance and depth grid points
#     tt: 2D array of values (e.g. travel times) at these points
#     shallowmode: Sets the extrapolation condition at shallow depth
#
#  < Returns
#     TTinterp: interpolation objection, used like tt = TTinterp(sdist,qdep)

function make_smtrace_table(dists,depths,tt,shallowmode="throw")
    
    # define interpolant object, with extrapolation options
    if shallowmode == "reflect" 
        TTinterp = LinearInterpolation((dists,depths),tt,
            extrapolation_bc=((Line(),Line()),(Reflect(),Line())))
    elseif shallowmode == "flat"
        TTinterp = LinearInterpolation((dists,depths),tt,
            extrapolation_bc=((Line(),Line()),(Flat(),Line())))
    elseif shallowmode == "line"
        TTinterp = LinearInterpolation((dists,depths),tt,
            extrapolation_bc=((Line(),Line()),(Line(),Line())))
    else
        TTinterp = LinearInterpolation((dists,depths),tt,
            extrapolation_bc=((Line(),Line()),(Throw(),Line())))
    end
    return TTinterp
end