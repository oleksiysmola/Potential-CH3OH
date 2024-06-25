using DataFrames
using LinearAlgebra

molecule = "CH3OH"
include("$(molecule).jl")

zMatrixFile = "z-matrix-$(molecule).txt"
nuclei = String[]
bondedToNucleus = Int64[]
bondLengths = Vector{Float64}()
angleToNucleus = Int64[]
bondAngles = Vector{Float64}()
dihedralToNucleus = Int64[]
dihedralAngles = Vector{Float64}()
nonRigid = Int64[]
nuclearMasses = Vector{Float64}()
println("Reading z-matrix...")
for line in eachline(zMatrixFile)
    splitLine = split(line, r"\s+")
    push!(nuclei, splitLine[2])
    push!(bondedToNucleus, parse(Int64, splitLine[3]))
    push!(bondLengths, parse(Float64, splitLine[4]))
    push!(angleToNucleus, parse(Int64, splitLine[5]))
    push!(bondAngles, parse(Float64, splitLine[6]))
    push!(dihedralToNucleus, parse(Int64, splitLine[7]))
    push!(dihedralAngles, parse(Float64, splitLine[8]))
    push!(nonRigid, parse(Int64, splitLine[9]))
    push!(nuclearMasses, parse(Float64, splitLine[10]))
end
zMatrix = DataFrame(Nuclei = nuclei, Bond = bondedToNucleus, Angle = angleToNucleus,
    Dihedral = dihedralToNucleus, NonRigid = nonRigid, Mass = nuclearMasses,     
)
println(zMatrix)
println("Obtained z-matrix!")

println("Defining reference geometry...")
numberOfAtoms = size(zMatrix)[1]
numberOfVibrationalModes = numberOfAtoms*3 - 6
localModeCoordinates = zeros(Float64, numberOfVibrationalModes)
convertToRadians = 2*pi/360
numberOfBondLengths = count(zMatrix[:, :Bond] .> 0)
numberOfAngles = count(zMatrix[:, :Angle] .> 0)
numberOfDihedrals = count(zMatrix[:, :Dihedral] .> 0)
localModeCoordinates[1:numberOfBondLengths] = bondLengths[2:numberOfAtoms]
localModeCoordinates[numberOfBondLengths+1:numberOfBondLengths+numberOfAngles] = bondAngles[3:numberOfAtoms].*convertToRadians
localModeCoordinates[numberOfBondLengths+numberOfAngles+1:numberOfVibrationalModes] = dihedralAngles[4:numberOfAtoms].*convertToRadians
cartesianCoordinates = computeCartesianCoordinates(localModeCoordinates)
println("Done!")

println("Obtaining symmetry operations...")
localModeSymmetryOperations = generateSymmetryOperationsLocalModeRepresentation()
println("Done!")

println("Generating initial potential parameters...")
potentialParameters = generateInitialPotentialParameters()
println(potentialParameters)
# for (label, parameter) in potentialParameters
#     println("$label  $parameter")
# end
println("Done!")
