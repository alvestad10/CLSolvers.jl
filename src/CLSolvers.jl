module CLSolvers

using LinearAlgebra
using Random
using Distributions
#using CUDA
using LinearAlgebra
using SparseArrays
using ForwardDiff
using ShiftedArrays, LabelledArrays

using DocStringExtensions


include("model.jl")
include("scheme.jl")
include("contour.jl")
include("contourDiscretization.jl")
include("analysis.jl")
include("plot_observables.jl")
include("sol_helpers.jl")


##################
#   EXPORT
##################
#export Plots


end
