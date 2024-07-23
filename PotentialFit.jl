using DataFrames
using Statistics
using StatsBase
using LsqFit
using GLM
using MLJLinearModels
using LinearAlgebra

molecule::String = "CH3OH"

include("$(molecule).jl")
include("TikhonovRegularisation.jl")

hartreeToWavenumberConversion::Float64 = 219474.63

collectedCoefficients::Array{Any} = []
if lowercase(readline()) == "parameters"
    dihedralMixLine::String = readline()
    splitDihedralMixLine::Array{SubString{String}} = split(dihedralMixLine, r"\s+")
    sameMultiMode::Vector{Vector{Int64}} = []
    if splitDihedralMixLine[end] != 0
        for i in 2:size(splitDihedralMixLine)[1]
            push!(sameMultiMode, parse.(Int64, split(String.(splitDihedralMixLine[i]), ",")))
        end
    end
    readingCoefficients::Bool = true
    while readingCoefficients
        coefficientInput::String = readline()
        if lowercase(coefficientInput) == "end"
            global readingCoefficients = false
            break
        end
        splitCoefficientInput::Array{SubString{String}} = split(coefficientInput, r"\s+")
        powers::Vector{Int64} = parse.(Int64, splitCoefficientInput[2:end - 1])
        multiModeVector::Vector{Bool} = powers[1:end-1] .> 0
        multiMode::Int64 = sum(multiModeVector)
        for multiModeGroup in sameMultiMode
            for i in 1:size(multiModeGroup)[1]
                multiModeOfGroup::Int64 = 0
                multiModeGroupVector::Vector{Int64} = zeros(size(multiModeGroup)[1])
                if multiModeVector[multiModeGroup[i]] > 0
                    multiModeOfGroup = size(multiModeGroup)[1]
                    multiModeGroupVector[i] = 1
                end
                multiMode = multiMode - sum(multiModeGroupVector) + multiModeOfGroup
            end
        end
        push!(collectedCoefficients, (String(splitCoefficientInput[1]), powers, parse(Float64, splitCoefficientInput[end]), Int64(sum(powers) - powers[end]), multiMode))
    end
end
expansionCoefficients::DataFrame = DataFrame(collectedCoefficients, [:torsionType, :expansionPowers, :coefficientValue, :expansionOrder, :multiMode])
expansionCoefficients = sort(expansionCoefficients, [:multiMode, :expansionOrder])

while lowercase(readline()) != "grid"
end

readingGrid::Bool = true
collectedGrid::Vector{Tuple} = []
while readingGrid
    gridLine::String = readline()
    if lowercase(gridLine) == "end"
        global readingGrid = false
        break
    end
    gridLineSplit::Array{SubString{String}} = split(gridLine, r"\s+")
    localModes::Vector{Float64} = parse.(Float64, gridLineSplit[1:end-2])
    push!(collectedGrid, (localModes, parse(Float64, gridLineSplit[end-1])))
end

grid::DataFrame = DataFrame(collectedGrid, [:localModes, :E])
numberOfModes::Int64 = size(grid[1, :localModes])[1]
numberOfPoints::Int64 = size(grid)[1]
println("Minimum in the potential energy:")
println(minimum(grid[:, :E]))
grid[:, :E] .= grid[:, :E].*hartreeToWavenumberConversion
grid[:, :E] .= grid[:, :E].-minimum(grid[:, :E])

# Weight factor by Partridge and Schwenke
function computeWeightOfPoint(energy::Float64, energyThreshold=15000.0::Float64)::Float64
    weight::Float64 = (tanh(−0.0006*(energy - energyThreshold)) + 1.002002002)/2.002002002
    if energy > 10000.0
        weight = weight/(0.0001*energy)
    else
        weight = weight/(0.0001*10000.0)
    end
    return weight
end

grid[:, :weight] .= computeWeightOfPoint.(grid[:, :E])
grid[:, :xi] .= generateXiCoordinates.(grid[:, :localModes])


expansionCoefficients = expansionCoefficients[expansionCoefficients.expansionOrder .< 4, :]
numberOfParameters::Int64 = size(expansionCoefficients)[1]
symmetryOperations::Array{Float64} = defineSymmetryOperations()
println("Invariance of current potential:")
xiCoordinates = grid[1, 1]
println(computePotentialEnergy(xiCoordinates, expansionCoefficients ,symmetryOperations))
xiCoordinates[1:end-1] = symmetryOperations[2, :, :]*xiCoordinates[1:end-1]
xiCoordinates[end] = xiCoordinates[end] + 2*pi/3
println(computePotentialEnergy(xiCoordinates, expansionCoefficients ,symmetryOperations))
xiCoordinates = grid[1, 1]
xiCoordinates[1:end-1] = symmetryOperations[3, :, :]*xiCoordinates[1:end-1]
xiCoordinates[end] = xiCoordinates[end] - 2*pi/3
println(computePotentialEnergy(xiCoordinates, expansionCoefficients ,symmetryOperations))
println("Initializing fit coordinates...")
@time xiPowers::Matrix{Float64} = setupFitVariables(grid, symmetryOperations, expansionCoefficients)
# @time grid::DataFrame = setupFitVariables(grid, symmetryOperations, expansionCoefficients)
println("Done!")


# Ridge Regression 
# potentialEnergyModel = @formula(E ~ .)
# reducedGrid::DataFrame = hcat(grid[:, :E], grid[:, 5:end])
# normalisedGrid::DataFrame
# for i in 1:numberOfParameters + 1
#     normalisedGrid[:, i] .= (reducedGrid[:, i] .- mean(reducedGrid[:, i]))./std(reducedGrid[:, i])
# end 
# lambda::Float64 = 1.0
# println("Begin fitting...")
# @time potentialFit = lm(potentialEnergyModel, reducedGrid, RidgeReg(lambda = lambda))
# println("Done!")

# fittedParameters::Vector{Float64} = coef(potentialFit)
# error = coeftable(potentialFit).se
# fittedPotential::Vector{Float64} = predict(potentialEnergyModel, reducedGrid)
# grid[:, :fittedEnergies] .= fittedPotential
# grid[:, :obsMinusCalc] .= grid[:, :E] .- grid[:, :fittedEnergies]

# Least Squares Fit
energies::Vector{Float64} = convert(Vector, grid[:, :E])
expansionParameters::Vector{Float64} = convert(Vector, expansionCoefficients[:, :coefficientValue])
weights::Vector{Float64} = convert(Vector, grid[:, :weight])
function potentialEnergyModel(xiPowers::Matrix{Float64}, expansionParameters::Vector{Float64})::Vector{Float64}
    numberOfPoints::Int64 = size(xiPowers)[1]
    potential::Vector{Float64} = zeros(numberOfPoints)
    for i in 1:numberOfPoints
        potential[i] = dot(xiPowers[i, :], expansionParameters)
    end
    return potential
end

function derivatives(xiPowers::Matrix{Float64}, expansionParameters::Vector{Float64})::Matrix{Float64}
    return xiPowers
end

normalisedEnergies::Vector{Float64} = (energies .- mean(energies))./std(energies)
normalisedXiPowers::Matrix{Float64} = zeros(numberOfPoints, numberOfParameters)
for i in 1:numberOfParameters
    normalisedXiPowers[:, i] = (xiPowers[:, i] .- mean(xiPowers[:, i]))./std(xiPowers[:, i]) 
end
println("Begin fitting...")
@time potentialFit = curve_fit((xiPowers, expansionParameters) -> potentialEnergyModel(xiPowers, expansionParameters),
    (xiPowers, expansionParameters) -> derivatives(xiPowers, expansionParameters),
    xiPowers, energies, weights, expansionParameters)
# @time fittedParameters::Vector{Float64} = MLJLinearModels.fit(RidgeRegression(), xiPowers[:, 2:end], energies)
println("Done!")

# println(potentialFit)
# println(size(expansionParameters))
# println(size(fittedParameters))

fittedParameters::Vector{Float64} = potentialFit.param
fittedPotential::Vector{Float64} = potentialEnergyModel(xiPowers, fittedParameters)
# fittedPotential::Vector{Float64} = potentialEnergyModel(xiPowers[:, 2:end], fittedParameters[2:end]) .+ fittedParameters[1]
grid[:, :fittedEnergies] .= fittedPotential
grid[:, :obsMinusCalc] .= grid[:, :E] .- grid[:, :fittedEnergies]
println("Parameters of fit:")
displayResult::String = ""
for i in 1:size(fittedParameters)[1]
    global displayResult = ""
    displayResult = displayResult*"$(expansionCoefficients[i, 1])"
    displayResult = displayResult*"  "
    powers = expansionCoefficients[i, 2]
    for j in 1:numberOfModes
        displayResult = displayResult*"$(powers[j])"*" "
    end
    displayResult = displayResult*"  "
    # displayResult = displayResult*"$(expansionCoefficients[i, 3])"*"  "
    displayResult = displayResult*"$(fittedParameters[i])"
    # displayResult = displayResult*"$(sigma[i])"
    println(displayResult)
end
println("Displaying energies at grid")
for i in 1:size(grid)[1]
    global displayResult = ""
    gridCoordinates = grid[i, 1]
    for j in 1:numberOfModes
        displayResult = displayResult*"$(gridCoordinates[j])"*" "
    end
    displayResult = displayResult*" "
    displayResult = displayResult*"$(grid[i, :E])"
    displayResult = displayResult*"  "
    displayResult = displayResult*"$(grid[i, :fittedEnergies])"
    displayResult = displayResult*"  "
    displayResult = displayResult*"$(grid[i, :obsMinusCalc])"
    println(displayResult)
end
println("rms:")
println(sqrt(mean(grid[:, :obsMinusCalc].^2)))
println("Error in parameters: ")
println(error)
