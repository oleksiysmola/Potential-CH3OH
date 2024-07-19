using DataFrames
using Statistics
using LinearAlgebra

molecule::String = "CH3OH"

include("$(molecule).jl")

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

expansionCoefficients = expansionCoefficients[expansionCoefficients.expansionOrder .< 3, :]
println(size(expansionCoefficients))
symmetryOperations::Array{Float64} = defineSymmetryOperations()
grid = setupFitVariables(grid, symmetryOperations, expansionCoefficients)
println(grid)
