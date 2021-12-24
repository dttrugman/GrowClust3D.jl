__precompile__()
module GrowClust3D

# external packages
using Printf, Random, Dates, DataFrames
using SharedArrays, Distributed
using Proj4: Transformation
using StatsBase: mean, median, mad, minimum, maximum
using Interpolations: LinearInterpolation, Flat, Line, Reflect
using DelimitedFiles: readdlm

# exports from these packages
export mean, median, mad, minimum, maximum

# exports from seismotrace
include("seismotrace.jl")
export read_vzmodel, find_moho, interp_vzmodel
export eflatten, eunflatten, trace_rays, first_arrivals
export write_table, smtrace_table 

# exports from inputs
include("inputs.jl")
export read_evlist, read_stlist, read_vzmodel, read_xcordata_proj
export check_gcinp, check_auxparams, read_gcinp

# exports from relocation
include("relocation.jl")
export evalrms, difclust1, difclust2, difclust3, map_distance, robomean, clustertree, xydist

end # module
