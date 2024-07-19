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
    return symmetryOperations
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
    
    xi::Vector{Float64} = zeros(Float64, 13)
    
    # Stretches
    xi[1] = 1.00 - exp(-a1*(localModeCoordinates[1] - rCOeq))
    xi[2] = 1.00 - exp(-a2*(localModeCoordinates[2] - rOHeq))
    xi[3] = 1.00 - exp(-b1*(localModeCoordinates[3] - rH1eq))
    xi[4] = 1.00 - exp(-b1*(localModeCoordinates[4] - rH1eq))
    xi[5] = 1.00 - exp(-b1*(localModeCoordinates[5] - rH1eq))

    # Angles
    xi[6] = cos(localModeCoordinates[6]) - cos(alphaCOHeq)
    xi[7] = cos(localModeCoordinates[7]) - cos(alphaOCHeq)
    xi[8] = cos(localModeCoordinates[8]) - cos(alphaOCHeq)
    xi[9] = cos(localModeCoordinates[9]) - cos(alphaOCHeq)

    # Define symmeterised dihedrals
    d32::Float64 = localModeCoordinates[11] - localModeCoordinates[12]
    d21::Float64 = localModeCoordinates[10] - localModeCoordinates[11]
    d13::Float64 = localModeCoordinates[12] - localModeCoordinates[10]
    xi[10] = (2*d32 - d21 - d13)/sqrt(6)
    xi[11] = (d21 - d13)/sqrt(2)
    
    # Torsion angle
    tau::Float64 = (localModeCoordinates[10] + localModeCoordinates[11] + localModeCoordinates[12])/3
    xi[12] = tau

    return xi
end

function setupFitVariables(grid::DataFrame, symmetryOperations::Array{Float64}, expansionCoefficients::DataFrame)::DataFrame
    numberOfSymmetryOperations::Int64 = size(symmetryOperations)[1]
    numberOfPoints::Int64 = size(grid)[1]
    # grid[:, :Xi] .= grid[:, :xi]
    # for i in 2:numberOfSymmetryOperations
    #     grid[:, :Xi] .= grid[:, :Xi] + symmetryOperations[i].*grid[:, :xi]
    # end
    # grid[:, :Xi] .= grid[:, :Xi]./6
    return grid
end

function computePotentialEnergy(xiCoordinates::Vector{Float64}, symmetryOperations::Array{Float64})::Float64
    potential = 0
    return potential
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