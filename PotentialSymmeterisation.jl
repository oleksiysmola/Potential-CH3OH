using DataFrames
using LinearAlgebra
using Symbolics
using Latexify

molecule = "CH3OH"
include("$(molecule).jl")

zMatrixFile = "z-matrix-$(molecule).txt"
nuclei = String[]
bondedToNucleus = Int64[]
bondLengths = String[]
angleToNucleus = Int64[]
bondAngles = String[]
dihedralToNucleus = Int64[]
dihedralAngles = String[]
nonRigid = Int64[]
nuclearMasses = Vector{Float64}()
println("Reading z-matrix...")
for line in eachline(zMatrixFile)
    splitLine = split(line, r"\s+")
    push!(nuclei, splitLine[2])
    push!(bondedToNucleus, parse(Int64, splitLine[3]))
    push!(bondLengths, splitLine[4])
    push!(angleToNucleus, parse(Int64, splitLine[5]))
    push!(bondAngles, splitLine[6])
    push!(dihedralToNucleus, parse(Int64, splitLine[7]))
    push!(dihedralAngles, splitLine[8])
    push!(nonRigid, parse(Int64, splitLine[9]))
    push!(nuclearMasses, parse(Float64, splitLine[10]))
end
zMatrix = DataFrame(Nuclei = nuclei, Bond = bondedToNucleus, Angle = angleToNucleus,
    Dihedral = dihedralToNucleus, NonRigid = nonRigid, Mass = nuclearMasses,     
)
println(zMatrix)
println("Obtained z-matrix!")
println("")
println("Defining local mode coordinates...")
numberOfAtoms = size(zMatrix)[1]
numberOfVibrationalModes = numberOfAtoms*3 - 6
convertToRadians = 2*pi/360
localModeCoordinates = Vector{Symbolics.Num}()
bondLengths = Symbol.(bondLengths[2:numberOfAtoms])
bondAngles = Symbol.(bondAngles[3:numberOfAtoms])
dihedralAngles = Symbol.(dihedralAngles[4:numberOfAtoms])
for bondLength in bondLengths
    @eval @variables $bondLength
    push!(localModeCoordinates, eval(bondLength))
end
for bondAngle in bondAngles
    @eval @variables $bondAngle
    push!(localModeCoordinates, eval(bondAngle))
end
for dihedral in dihedralAngles
    @eval @variables $dihedral
    push!(localModeCoordinates, eval(dihedral))
end
println(localModeCoordinates)
println(localModeCoordinates[3] + localModeCoordinates[4])
println("Done!")
println("")

println("Obtaining symmetry operations...")
localModeSymmetryOperations = generateSymmetryOperationsLocalModeRepresentation()
println("Done!")
println("")

println("Generating initial potential parameters...")
potentialParameters = generateInitialPotentialParameters()
println(potentialParameters)
# for (label, parameter) in potentialParameters
#     println("$label  $parameter")
# end
println("Done!")
x = Symbol("r_1")
@eval @variables $x
println(r_1)
# println(expr2)

