# Test 1D problem with ray tracing
println("\nRunning ray-tracing test case (1D velocity model).")
run(`julia -t1 ../examples/run_growclust3D.jl test.trace1D.inp`)
println("\nCOMPLETED: ray-tracing test case (1D velocity model).")
println()

# Test 1D problem with nonlinloc grids
println("\nRunning nonlinloc grid test case (1D velocity model).")
run(`julia -t1 ../examples/run_growclust3D.jl test.grid1D.inp`)
println("\nCOMPLETED: nonlinloc grid test case (1D velocity model).")
println()

# Test 3D problem with nonlinloc grids
println("\nRunning nonlinloc grid test case (3D velocity model).")
run(`julia -t1 ../examples/run_growclust3D.jl test.grid3D.inp`)
println("\nCOMPLETED: nonlinloc grid test case (3D velocity model).")
println()