function LsqCurveFit(potentialEnergyModel::Function, derivatives::Function, xiPowers::Matrix{Float64}, expansionParameters::Vector{Float64}, energies::Vector{Float64}, weights::Vector{Float64})::Vector{Float64}
    @time potentialFit = curve_fit((xiPowers, expansionParameters) -> potentialEnergyModel(xiPowers, expansionParameters),
    (xiPowers, expansionParameters) -> derivatives(xiPowers, expansionParameters),
    xiPowers, energies, weights, expansionParameters)
    return potentialFit.param
end

function TikhonovRegularisation(x::Matrix{Float64}, parameters::Vector{Float64}, y::Vector{Float64}, regularisationParameter=1.0::Float64)::Vector{Float64}
    numberOfParameters::Int64 = size(parameters)[1]
    tikhonovMatrix::Matrix{Float64} = Matrix(1I, numberOfParameters, numberOfParameters)
    solvedCoefficients::Vector{Float64} = inv((transpose(x)*x - transpose(tikhonovMatrix)*tikhonovMatrix))*transpose(x)*y
    return solvedCoefficients
end