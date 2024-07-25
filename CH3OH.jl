function defineSymmetryOperations(case="C3v(M)"::String)::Array{Float64}
    if case == "C3v(M)"
        symmetryOperations::Array{Float64} = zeros(Int64, 6, 11, 11)
        symmetryOperations[1, :, :] = Matrix(1I, 11, 11)
        # (123)
        symmetryOperations[2, :, :] = [
            1 0 0 0 0 0 0 0 0 0 0;
            0 1 0 0 0 0 0 0 0 0 0;
            0 0 0 1 0 0 0 0 0 0 0; 
            0 0 0 0 1 0 0 0 0 0 0; 
            0 0 1 0 0 0 0 0 0 0 0;
            0 0 0 0 0 1 0 0 0 0 0;
            0 0 0 0 0 0 0 1 0 0 0;
            0 0 0 0 0 0 0 0 1 0 0;
            0 0 0 0 0 0 1 0 0 0 0;
            0 0 0 0 0 0 0 0 0 -1/2 sqrt(3)/2;
            0 0 0 0 0 0 0 0 0 -sqrt(3)/2 -1/2;            
            ]
        # (132)
        symmetryOperations[3, :, :] = [
            1 0 0 0 0 0 0 0 0 0 0;
            0 1 0 0 0 0 0 0 0 0 0;
            0 0 0 0 1 0 0 0 0 0 0;
            0 0 1 0 0 0 0 0 0 0 0;
            0 0 0 1 0 0 0 0 0 0 0;
            0 0 0 0 0 1 0 0 0 0 0;
            0 0 0 0 0 0 0 0 1 0 0;
            0 0 0 0 0 0 1 0 0 0 0;
            0 0 0 0 0 0 0 1 0 0 0;
            0 0 0 0 0 0 0 0 0 -1/2 -sqrt(3)/2;
            0 0 0 0 0 0 0 0 0 sqrt(3)/2 -1/2;
        ]
        # (12)*
        symmetryOperations[4, :, :] = [
            1 0 0 0 0 0 0 0 0 0 0;
            0 1 0 0 0 0 0 0 0 0 0;
            0 0 0 1 0 0 0 0 0 0 0;
            0 0 1 0 0 0 0 0 0 0 0;
            0 0 0 0 1 0 0 0 0 0 0;
            0 0 0 0 0 1 0 0 0 0 0;
            0 0 0 0 0 0 0 1 0 0 0;
            0 0 0 0 0 0 1 0 0 0 0;
            0 0 0 0 0 0 0 0 1 0 0;
            0 0 0 0 0 0 0 0 0 -1/2 sqrt(3)/2;
            0 0 0 0 0 0 0 0 0 sqrt(3)/2 1/2;
        ]
        # (23)*
        symmetryOperations[5, :, :] = [
            1 0 0 0 0 0 0 0 0 0 0;
            0 1 0 0 0 0 0 0 0 0 0;    
            0 0 1 0 0 0 0 0 0 0 0;
            0 0 0 0 1 0 0 0 0 0 0;
            0 0 0 1 0 0 0 0 0 0 0;
            0 0 0 0 0 1 0 0 0 0 0;
            0 0 0 0 0 0 1 0 0 0 0;
            0 0 0 0 0 0 0 0 1 0 0;
            0 0 0 0 0 0 0 1 0 0 0;
            0 0 0 0 0 0 0 0 0 1 0;
            0 0 0 0 0 0 0 0 0 0 -1;
        ]
        # (13)*
        symmetryOperations[6, :, :] = [
            1 0 0 0 0 0 0 0 0 0 0;
            0 1 0 0 0 0 0 0 0 0 0;    
            0 0 0 0 1 0 0 0 0 0 0;
            0 0 0 1 0 0 0 0 0 0 0;
            0 0 1 0 0 0 0 0 0 0 0;
            0 0 0 0 0 1 0 0 0 0 0;
            0 0 0 0 0 0 0 0 1 0 0;
            0 0 0 0 0 0 0 1 0 0 0;
            0 0 0 0 0 0 1 0 0 0 0;
            0 0 0 0 0 0 0 0 0 -1/2 -sqrt(3)/2;
            0 0 0 0 0 0 0 0 0 -sqrt(3)/2 1/2;
        ]
    end
    return symmetryOperations
end

function defineTorsionSymmetryOperations(fourierOrder::Int64, case="C3v(M)"::String)::Array{Float64}
    if case == "C3v(M)"
        torsionSymmetryOperations::Array{Float64} = zeros(Int64, 6, 2, 2)
        torsionSymmetryOperations[1, :, :] = Matrix(1I, 2, 2)
        # (123)
        torsionSymmetryOperations[2, :, :] = [
            cos(2*pi*fourierOrder/3) -sin(2*pi*fourierOrder/3);
            sin(2*pi*fourierOrder/3) cos(2*pi*fourierOrder/3);           
            ]
        # (132)
        torsionSymmetryOperations[3, :, :] = [
            cos(2*pi*fourierOrder/3) sin(2*pi*fourierOrder/3);
            -sin(2*pi*fourierOrder/3) cos(2*pi*fourierOrder/3);           
            ]
        # (12)*
        torsionSymmetryOperations[4, :, :] = [
            cos(2*pi*fourierOrder/3) -sin(2*pi*fourierOrder/3);
            -sin(2*pi*fourierOrder/3) -cos(2*pi*fourierOrder/3);           
            ]
        # (23)*
        torsionSymmetryOperations[5, :, :] = [
            1 0;
            0 -1;           
            ]
        # (13)*
        torsionSymmetryOperations[6, :, :] = [
            cos(2*pi*fourierOrder/3) sin(2*pi*fourierOrder/3);
            sin(2*pi*fourierOrder/3) -cos(2*pi*fourierOrder/3);           
            ]
    end
    return torsionSymmetryOperations
end

function generateXiCoordinates(localModeCoordinates::Vector{Float64})::Vector{Float64}
    convertToRadians::Float64 = 2*pi/360
    # Define equilibrium parameters
    rCOeq::Float64 = 1.4296
    rOHeq::Float64 = 0.95887
    rH1eq::Float64 = 1.092294
    alphaCOHeq::Float64 = 107.9812*convertToRadians
    alphaOCHeq::Float64 = 110.6646*convertToRadians
    a1::Float64 = 1.44
    a2::Float64 = 1.66
    b1::Float64 = 1.55
    xi::Vector{Float64} = zeros(Float64, 12)
    
    # Stretches
    xi[1] = 1.00 - exp(-a1*(localModeCoordinates[1] - rCOeq))
    xi[2] = 1.00 - exp(-a2*(localModeCoordinates[2] - rOHeq))
    xi[3] = 1.00 - exp(-b1*(localModeCoordinates[3] - rH1eq))
    xi[4] = 1.00 - exp(-b1*(localModeCoordinates[4] - rH1eq))
    xi[5] = 1.00 - exp(-b1*(localModeCoordinates[5] - rH1eq))

    # Angles
    xi[6] = cos(localModeCoordinates[6]*convertToRadians) - cos(alphaCOHeq)
    xi[7] = localModeCoordinates[7]*convertToRadians - alphaOCHeq
    xi[8] = localModeCoordinates[8]*convertToRadians - alphaOCHeq
    xi[9] = localModeCoordinates[9]*convertToRadians - alphaOCHeq

    # Define symmeterised dihedrals
    d12eq::Float64 = 2.165647870654562
    d23eq::Float64 = 1.951889565870463
    d13eq::Float64 = 2.165647870654562
    d12::Float64 = mod((localModeCoordinates[11] - localModeCoordinates[10])*convertToRadians + 2*pi, 2*pi)
    d23::Float64 = mod((localModeCoordinates[12] - localModeCoordinates[11])*convertToRadians + 2*pi, 2*pi)
    d13::Float64 = mod((localModeCoordinates[10] - localModeCoordinates[12])*convertToRadians + 2*pi, 2*pi)
    xi[10] = (2*d23 - d13 - d12)/sqrt(6)
    xi[11] = (d13 - d12)/sqrt(2)
    
    # Torsion angle
    tau::Float64 = ((localModeCoordinates[10] + localModeCoordinates[11] + localModeCoordinates[12])*convertToRadians - 2*pi)/3
    xi[12] = tau

    return xi
end

# function setupFitVariables(grid::DataFrame, symmetryOperations::Array{Float64}, expansionCoefficients::DataFrame)::DataFrame
#     numberOfSymmetryOperations::Int64 = size(symmetryOperations)[1]
#     numberOfPoints::Int64 = size(grid)[1]
#     numberOfExpansionCoefficients::Int64 = size(expansionCoefficients)[1]
#     xiPowers::Matrix{Float64} = zeros(Float64, numberOfPoints, numberOfExpansionCoefficients)
#     for i in 1:numberOfPoints
#         xi::Vector{Float64} = grid[i, :xi]
#         for j in 1:numberOfExpansionCoefficients
#             torsionType::String = expansionCoefficients[j, :torsionType]
#             powers::Vector{Int64} = expansionCoefficients[j, :expansionPowers]
#             torsionSymmetryOperations::Array{Float64} = defineTorsionSymmetryOperations(powers[end])
#             torsionMode::Vector{Float64} = [cos(powers[end]*xi[end]), sin(powers[end]*xi[end])]
#             for k in 1:numberOfSymmetryOperations
#                 transformedTorsionMode::Vector{Float64} = torsionSymmetryOperations[k, :, :]*torsionMode
#                 transformedXi::Vector{Float64} = symmetryOperations[k, :, :]*xi[1:end-1]
#                 if torsionType == "f"
#                     xiPowers[i, j] = xiPowers[i, j] + transformedTorsionMode[1]*prod(transformedXi.^powers[1:end-1])
#                 else
#                     xiPowers[i, j] = xiPowers[i, j] + transformedTorsionMode[2]*prod(transformedXi.^powers[1:end-1])
#                 end
#             end
#         end
#     end
#     grid = hcat(grid, DataFrame(xiPowers, :auto))
#     return grid
# end

function setupFitVariables(grid::DataFrame, symmetryOperations::Array{Float64}, expansionCoefficients::DataFrame)::Matrix{Float64}
    numberOfSymmetryOperations::Int64 = size(symmetryOperations)[1]
    numberOfPoints::Int64 = size(grid)[1]
    numberOfExpansionCoefficients::Int64 = size(expansionCoefficients)[1]
    xiPowers::Matrix{Float64} = zeros(Float64, numberOfPoints, numberOfExpansionCoefficients)
    for i in 1:numberOfPoints
        xi::Vector{Float64} = grid[i, :xi]
        for j in 1:numberOfExpansionCoefficients
            torsionType::String = expansionCoefficients[j, :torsionType]
            powers::Vector{Int64} = expansionCoefficients[j, :expansionPowers]
            torsionSymmetryOperations::Array{Float64} = defineTorsionSymmetryOperations(powers[end])
            torsionMode::Vector{Float64} = [cos(powers[end]*xi[end]), sin(powers[end]*xi[end])]
            for k in 1:numberOfSymmetryOperations
                transformedTorsionMode::Vector{Float64} = torsionSymmetryOperations[k, :, :]*torsionMode
                transformedXi::Vector{Float64} = symmetryOperations[k, :, :]*xi[1:end-1]
                if torsionType == "f"
                    xiPowers[i, j] = xiPowers[i, j] + transformedTorsionMode[1]*prod(transformedXi.^powers[1:end-1])
                else
                    xiPowers[i, j] = xiPowers[i, j] + transformedTorsionMode[2]*prod(transformedXi.^powers[1:end-1])
                end
            end
        end
    end
    return xiPowers
end

function computePotentialEnergy(xiCoordinates::Vector{Float64}, expansionCoefficients::DataFrame, symmetryOperations::Array{Float64})::Float64
    potential::Float64 = 0.0
    numberOfParameters::Int64 = size(expansionCoefficients)[1]
    numberOfSymmetryOperations::Int64 = size(symmetryOperations)[1]
    xiPowers::Vector{Float64} = zeros(numberOfParameters)
    for i in 1:numberOfParameters
        torsionType::String = expansionCoefficients[i, 1]
        powers::Vector{Int64} = expansionCoefficients[i, 2]
        torsionSymmetryOperations::Array{Float64} = defineTorsionSymmetryOperations(powers[end])
        torsionMode::Vector{Float64} = [cos(powers[end]*xiCoordinates[end]), sin(powers[end]*xiCoordinates[end])]
        for j in 1:numberOfSymmetryOperations
            transformedTorsionMode::Vector{Float64} = torsionSymmetryOperations[j, :, :]*torsionMode
            transformedXi::Vector{Float64} = symmetryOperations[j, :, :]*xiCoordinates[1:end-1]
            if torsionType == "f"
                xiPowers[i] = xiPowers[i] + transformedTorsionMode[1]*prod(transformedXi.^powers[1:end-1])
            else
                xiPowers[i] = xiPowers[i] + transformedTorsionMode[2]*prod(transformedXi.^powers[1:end-1])
            end
        end
    end
    potential = dot(xiPowers, expansionCoefficients[:, 3])
    return potential
end

function checkPotentialForInvariance(grid::DataFrame, expansionCoefficients::DataFrame, symmetryOperations::Array{Float64})
    xiCoordinates::Vector{Float64} = copy(grid[1, 1])
    potentialBeforeTransformation::Float64 = computePotentialEnergy(xiCoordinates, expansionCoefficients, symmetryOperations)
    xiCoordinates[1:end-1] = symmetryOperations[1, :, :]*xiCoordinates[1:end-1]
    println("Symmetry Operation E: ", computePotentialEnergy(xiCoordinates, expansionCoefficients, symmetryOperations) - potentialBeforeTransformation)
    xiCoordinates[1:end-1] = symmetryOperations[2, :, :]*xiCoordinates[1:end-1]
    xiCoordinates[end] = xiCoordinates[end] + 2*pi/3
    println("Symmetry Operation (123): ", computePotentialEnergy(xiCoordinates, expansionCoefficients, symmetryOperations) - potentialBeforeTransformation)
    xiCoordinates = copy(grid[1, 1])
    xiCoordinates[1:end-1] = symmetryOperations[3, :, :]*xiCoordinates[1:end-1]
    xiCoordinates[end] = xiCoordinates[end] - 2*pi/3
    println("Symmetry Operation (132): ", computePotentialEnergy(xiCoordinates, expansionCoefficients, symmetryOperations) - potentialBeforeTransformation)
    xiCoordinates = copy(grid[1, 1])
    xiCoordinates[1:end-1] = symmetryOperations[4, :, :]*xiCoordinates[1:end-1]
    xiCoordinates[end] = -xiCoordinates[end] - 2*pi/3
    println("Symmetry Operation (12)*: ", computePotentialEnergy(xiCoordinates, expansionCoefficients, symmetryOperations) - potentialBeforeTransformation)
    xiCoordinates = copy(grid[1, 1])
    xiCoordinates[1:end-1] = symmetryOperations[5, :, :]*xiCoordinates[1:end-1]
    xiCoordinates[end] = -xiCoordinates[end]
    println("Symmetry Operation (23)*: ", computePotentialEnergy(xiCoordinates, expansionCoefficients, symmetryOperations) - potentialBeforeTransformation)
    xiCoordinates = copy(grid[1, 1])
    xiCoordinates[1:end-1] = symmetryOperations[6, :, :]*xiCoordinates[1:end-1]
    xiCoordinates[end] = -xiCoordinates[end] + 2*pi/3
    println("Symmetry Operation (13)*: ", computePotentialEnergy(xiCoordinates, expansionCoefficients, symmetryOperations) - potentialBeforeTransformation)
end

function computeCartesianCoordinates(localModeCoordinates::Vector{Float64}, case="bond-fixed-dihedral"::String)::Matrix{Float64}
    cartesianCoordinates::Matrix{Float64} = zeros(6, 3)
    if case == "bond-fixed-dihedral"
        # O
        cartesianCoordinates[2, 3] = localModeCoordinates[1]
        # H
        cartesianCoordinates[3, 1] = sin(localModeCoordinates[6])*localModeCoordinates[2]
        cartesianCoordinates[3, 3] = cartesianCoordinates[2, 3] - cos(localModeCoordinates[6])*localModeCoordinates[2]
        # H1
        cartesianCoordinates[4, 1] = sin(localModeCoordinates[7])*cos(localModeCoordinates[10])*localModeCoordinates[3]
        cartesianCoordinates[4, 2] = sin(localModeCoordinates[7])*sin(localModeCoordinates[10])*localModeCoordinates[3]
        cartesianCoordinates[4, 3] = cos(localModeCoordinates[7])*localModeCoordinates[3]
        # H2
        cartesianCoordinates[5, 1] = sin(localModeCoordinates[8])*cos(localModeCoordinates[11])*localModeCoordinates[4]
        cartesianCoordinates[5, 2] = sin(localModeCoordinates[8])*sin(localModeCoordinates[11])*localModeCoordinates[4]
        cartesianCoordinates[5, 3] = cos(localModeCoordinates[8])*localModeCoordinates[4]
        # H3
        cartesianCoordinates[6, 1] = sin(localModeCoordinates[9])*cos(localModeCoordinates[12])*localModeCoordinates[5]
        cartesianCoordinates[6, 2] = sin(localModeCoordinates[9])*sin(localModeCoordinates[12])*localModeCoordinates[5]
        cartesianCoordinates[6, 3] = cos(localModeCoordinates[9])*localModeCoordinates[5]
    end
    return cartesianCoordinates
end