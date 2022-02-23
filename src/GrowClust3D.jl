__precompile__()
module GrowClust3D

# external packages
using Printf, Random, Dates, DataFrames
using SharedArrays, Distributed
using Proj4: Transformation
using StatsBase: mean, median, mad, minimum, maximum, percentile, sample, std
using Interpolations: LinearInterpolation, CubicSplineInterpolation
using Interpolations: Flat, Line, Reflect, Throw
using DelimitedFiles: readdlm

# exports from these packages
export mean, median, mad, minimum, maximum, percentile, sample, std

# exports from seismotrace
include("seismotrace.jl")
export read_vzmodel, find_moho, interp_vzmodel
export eflatten, eunflatten, trace_rays, first_arrivals
export write_smtrace_table, read_smtrace_table, make_smtrace_table

# exports from nllgrid
include("nllgrid.jl")
export read_nll_head, make_nll_interp, check_proj

# exports from inputs
include("inputs.jl")
export read_evlist, read_stlist, read_vzmodel, read_xcordata_proj
export check_gcinp, check_auxparams, read_gcinp

# exports from relocation
include("relocation1D.jl")
export evalrms, difclust1, difclust2, difclust3, map_distance, robomean, clustertree, xydist
include("relocation3D.jl")
export difclust1_3D, difclust2_3D, difclust3_3D, clustertree_3D

end # module
