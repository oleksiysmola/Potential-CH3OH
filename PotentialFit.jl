using DataFrames
using LinearAlgebra
using SymPy

molecule::String = "CH3OH"

include("$(molecule).jl")

println(defineSymmetryOperations()[2, :, :])
