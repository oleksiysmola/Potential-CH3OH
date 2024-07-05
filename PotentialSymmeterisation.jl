using DataFrames
using LinearAlgebra
using SymPy
using Distributed
using TensorOperations

# Access the dynamically created symbols
x = SymPy.Sym("r_1")
println(typeof(x))
println(SymPy.Sym("0") + SymPy.Sym("1"))

# Create and manipulate expressions using the dynamic symbols
expr = x^2  - 1 
# println(typeof(Eq(expr)))
println("Expression: ", expr)
println(diff(expr, x))
println(solve(expr, x))

molecule::String = "CH3OH"
include("$(molecule).jl")


# @everywhere using DataFrames
# @everywhere using LinearAlgebra
# @everywhere using SymPy


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
localModeCoordinates::Vector{SymPy.Sym} = Vector{SymPy.Sym}()
bondLengthSymbols::Vector{SymPy.Sym} = SymPy.Sym.(bondLengths[2:numberOfAtoms])
bondAnglesSymbols::Vector{SymPy.Sym} = SymPy.Sym.(bondAngles[3:numberOfAtoms])
dihedralAnglesSymbols::Vector{SymPy.Sym} = SymPy.Sym.(dihedralAngles[4:numberOfAtoms])
for bondLength in bondLengthSymbols
    push!(localModeCoordinates, bondLength)
end
for bondAngle in bondAnglesSymbols
    push!(localModeCoordinates, bondAngle)
end
for dihedral in dihedralAnglesSymbols
    push!(localModeCoordinates, dihedral)
end
expansionCoefficientsTensorForm, coefficients = setupTensorFormPotential()
println(expansionCoefficientsTensorForm)
println(coefficients)
println(localModeCoordinates)
println("Done!")
println("")

println("Obtaining symmetry operations...")
symmetryOperations::Array{SymPy.Sym} = generateSymmetryOperationsXi()
@time solveForSymmetricCoefficients(expansionCoefficientsTensorForm, coefficients, symmetryOperations)
println("Done!")
println("")

println("Generating xi coordinates...")
xi::Vector{SymPy.Sym} = generateXiCoordinates(localModeCoordinates)
println(xi)
println("Done!")
println("")

println("Generating initial potential parameters...")
@time potentialParameters::DataFrame = generateInitialPotentialParameters()
println(potentialParameters)
println("Done!")
println("")

println("Obtain all xi power terms...")
@time potentialTerms::DataFrame = obtainTransformedPotentialTermsLocalMode(potentialParameters, localModeCoordinates, symmetryOperations)
println(potentialTerms)
println("Done!")
println("")
# println("Computing potentials...")
# numberOfSymmetryOperations::Int64 = size(symmetryOperations)[1]
# potentials::Vector{SymPy.Sym} = zeros(numberOfSymmetryOperations + 1)
# for i in 1:size(potentialTerms)[1]
#     for j in 1:size(potentials)[1]
#         potentials[j] = potentials[j] + potentialTerms[i, names(potentialTerms)[1]]*potentialTerms[i, names(potentialTerms)[j + numberOfVibrationalModes + 2]]
#     end
# end
# println("Done!")
# println("")

# println("Taking difference between transformed potentials...")
# equations::Vector{SymPy.Sym} = Vector{SymPy.Sym}()
# for i in 2:size(potentials)[1]
#     push!(equations, potentials[i] - potentials[1])
# end
# println("Done!")
# println("")


# println("Solving for coefficients...")
# solutions = solve(equations, potentialTerms[:, :Labels])
# @everywhere function solveEquation(equation::SymPy.Sym, variables::Vector{SymPy.Sym})
#     return solve(equation, variables)
# end
# @time solutions = pmap(equation -> solveEquation(equation, potentialTerms[:, :Labels]), equations)
println(solutions)
println("Done!")
println("")