# Testing options
test_trace1D = true
test_grid1D = false # requires NLL grids
test_grid3D = false # requires NLL grids

# Test 1D problem with ray tracing
if test_trace1D
    println("\nRunning ray-tracing test case (1D velocity model).")
    run(`julia -t1 run_growclust3D.jl test.trace1D.inp`)
    println("\nCOMPLETED: ray-tracing test case (1D velocity model).")
    println()
end

# Test 1D problem with nonlinloc grids
if test_grid1D
    println("\nRunning nonlinloc grid test case (1D velocity model).")
    run(`julia -t1 run_growclust3D.jl test.grid1D.inp`)
    println("\nCOMPLETED: nonlinloc grid test case (1D velocity model).")
    println()
end

# Test 3D problem with nonlinloc grids
if test_grid3D
    println("\nRunning nonlinloc grid test case (3D velocity model).")
    run(`julia -t1 run_growclust3D.jl test.grid3D.inp`)
    println("\nCOMPLETED: nonlinloc grid test case (3D velocity model).")
    println()
end