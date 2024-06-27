using DataFrames
using LinearAlgebra
using Symbolics

molecule::String = "CH3OH"
include("$(molecule).jl")

zMatrixFile::String = "z-matrix-$(molecule).txt"
nuclei::Array{String} = String[]
bondedToNucleus::Array{Int64} = Int64[]
bondLengths::Array{String} = String[]
angleToNucleus::Array{Int64} = Int64[]
bondAngles::Array{String} = String[]
dihedralToNucleus::Array{Int64} = Int64[]
dihedralAngles::Array{String} = String[]
nonRigid::Array{Int64} = Int64[]
nuclearMasses::Array{Float64} = Vector{Float64}()
println("Reading z-matrix...")
for line in eachline(zMatrixFile)
    splitLine::Array{SubString{String}, 1} = split(line, r"\s+")
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
zMatrix::DataFrame = DataFrame(Nuclei = nuclei, Bond = bondedToNucleus, Angle = angleToNucleus,
    Dihedral = dihedralToNucleus, NonRigid = nonRigid, Mass = nuclearMasses,     
)
println(zMatrix)
println("Obtained z-matrix!")
println("")
println("Defining local mode coordinates...")
numberOfAtoms::Int64 = size(zMatrix)[1]
numberOfVibrationalModes::Int64 = numberOfAtoms*3 - 6
localModeCoordinates::Vector{Symbolics.Num} = Vector{Symbolics.Num}()
bondLengthSymbols::Vector{Symbol} = Symbol.(bondLengths[2:numberOfAtoms])
bondAnglesSymbols::Vector{Symbol} = Symbol.(bondAngles[3:numberOfAtoms])
dihedralAnglesSymbols::Vector{Symbol} = Symbol.(dihedralAngles[4:numberOfAtoms])
for bondLength in bondLengthSymbols
    @eval @variables $bondLength
    push!(localModeCoordinates, eval(bondLength))
end
for bondAngle in bondAnglesSymbols
    @eval @variables $bondAngle
    push!(localModeCoordinates, eval(bondAngle))
end
for dihedral in dihedralAnglesSymbols
    @eval @variables $dihedral
    push!(localModeCoordinates, eval(dihedral))
end
println("Done!")
println("")

println("Obtaining symmetry operations...")
symmetryOperations::Array{Int64} = generateSymmetryOperationsLocalModeRepresentation()
println("Done!")
println("")

println("Generating xi coordinates...")
xi::Vector{Symbolics.Num} = generateXiCoordinates(localModeCoordinates)
println(xi)
println("Done!")
println("")

println("Generating initial potential parameters...")
@time potentialParameters::DataFrame = generateInitialPotentialParameters()
println(potentialParameters)
println("Done!")
println("")

println("Symmetrerise potential...")
@time potentialTerms::DataFrame = obtainTransformedPotentialTermsLocalMode(potentialParameters, localModeCoordinates, symmetryOperations)
println(potentialTerms)
numberOfSymmetryOperations::Int64 = size(symmetryOperations)[1]
potentials::Vector{Symbolics.Num} = zeros(numberOfSymmetryOperations + 1)
for i in 1:size(potentialTerms)[1]
    for j in 1:size(potentials)[1]
        potentials[j] = potentials[j] + potentialTerms[i, names(potentialTerms)[1]]*potentialTerms[i, names(potentialTerms)[j + numberOfVibrationalModes + 2]]
    end
end
equations::Vector{Equation} = Vector{Equation}()
for i in 2:size(potentials)[1]
    push!(equations, potentials[i] - potentials[1] ~ 0)
end
@time solutions = Symbolics.solve_for(equations, potentialTerms[:, :Labels])
println(solutions)
println("Done!")
println("")
