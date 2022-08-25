### Script To Generate NonLinLoc Grid Files ####
# Daniel Trugman, 2022

### Modules ###

using Printf

### Setup ###

## NonLinLoc Source code (must be precompiled)
#  - assumes the NonLinLoc bin directory is on your path
#  - if not, just add the full path before Vel2Grid and Grid2Time
srcV2G = "Vel2Grid"
srcG2T = "Grid2Time"

## Directories
workDIR = "./data/nll/"
runDIR = workDIR * "nll-run/"
grdDIR = workDIR * "nll-grd/"
grdDIR3 = workDIR * "nll-grd3/"

## NonLinLoc Run Files (separated for clarity)
# 1D models
fV2G = runDIR * "v2g.inp"
fG2TP = runDIR * "g2t.P.inp"
fG2TS = runDIR * "g2t.S.inp"
# 3D models
fV2G3 = runDIR * "v2g3.inp" # 3D model
fG2TP3 = runDIR * "g2t3.P.inp" # 3D grid
fG2TS3 = runDIR * "g2t3.S.inp" # 3D grid

## Script options
make1D = true # option to make 1D travel time tables (quick, low memory requirement)
make3D = true # option to make 3D travel time grids (takes longer, requires more memory)
makeS = false # we actually don't need to make S-wave tables for constant Vp/Vs

### Generate 1D Tables

if make1D

    # Make output path
    println("\nWorking on 1D grids, output path: ", grdDIR)
    mkpath(grdDIR)

    # Vel2Grid
    println("\nRunning Vel2Grid...")
    run(`$srcV2G $fV2G`)

    # Grid2Time: P-wave
    println("\nRunning Grid2Time: P-wave")
    run(`$srcG2T $fG2TP`)

    # Grid2Time: S-wave
    if makeS # optional
        println("\nRunning Grid2Time: S-wave")
        run(`$srcG2T $fG2TS`)
    end

    # Done
    println("\nCompleted Grid Creation")

end

### Generate 1D Tables

if make3D

    # Make output path
    println("\nWorking on 3D grids, output path: ", grdDIR3)
    mkpath(grdDIR3)

    # Vel2Grid
    println("\nRunning Vel2Grid...")
    run(`$srcV2G $fV2G3`)

    # Grid2Time: P-wave
    println("\nRunning Grid2Time: P-wave")
    run(`$srcG2T $fG2TP3`)

    # Grid2Time: S-wave
    if makeS # optional
        println("\nRunning Grid2Time: S-wave")
        run(`$srcG2T $fG2TS3`)
    end

    # Done
    println("\nCompleted Grid Creation")

end

