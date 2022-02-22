# Test 1D problem with ray tracing
println("\nRunning ray-tracing test case (1D velocity model).")
run(`julia -t1 test_ssprings1D.jl test.trace.inp`)
println("\nCOMPLETED: ray-tracing test case (1D velocity model).")
println()

# Test 1D problem with nonlinloc grids
println("\nRunning nonlinloc grid test case (1D velocity model).")
run(`julia -t1 test_ssprings1D.jl test.grid.inp`)
println("\nCOMPLETED: nonlinloc grid test case (1D velocity model).")
println()