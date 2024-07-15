function computeCartesianCoordinates(localModeCoordinates::Vector{SymPy.Sym}, case="bond-fixed-dihedral"::String)
    cartesianCoordinates = zeros(6, 3)
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

function generateSymmetryOperationsLocalModeRepresentation(case="C3v(M)"::String)::Array{Int64}
    permuation123::Matrix{Int64} = Int64[0 1 0; 0 0 1; 1 0 0]
    permuation132::Matrix{Int64} = Int64[0 0 1; 1 0 0; 0 1 0]
    permutationInversion12::Matrix{Int64} = Int64[0 1; 1 0]
    permutationInversion23::Matrix{Int64} = Int64[0 1; 1 0]
    permutationInversion13::Matrix{Int64} = Int64[0 0 1; 0 1 0; 1 0 0]
    if case == "C3v(M)"
        localModeSymmetryOperations::Array{Int64} = zeros(Int64, 6, 12, 12)
        for i in 1:6
            localModeSymmetryOperations[i, :, :] = Matrix(1I, 12, 12)
        end
        localModeSymmetryOperations[2, 3:5, 3:5] = permuation123
        localModeSymmetryOperations[2, 7:9, 7:9] = permuation123
        localModeSymmetryOperations[2, 10:12, 10:12] = permuation123
        localModeSymmetryOperations[3, 3:5, 3:5] = permuation132
        localModeSymmetryOperations[3, 7:9, 7:9] = permuation132
        localModeSymmetryOperations[3, 10:12, 10:12] = permuation132
        localModeSymmetryOperations[4, 3:4, 3:4] = permutationInversion12
        localModeSymmetryOperations[4, 7:8, 7:8] = permutationInversion12
        localModeSymmetryOperations[4, 10:11, 10:11] = permutationInversion12
        localModeSymmetryOperations[5, 4:5, 4:5] = permutationInversion23
        localModeSymmetryOperations[5, 8:9, 8:9] = permutationInversion23
        localModeSymmetryOperations[5, 11:12, 11:12] = permutationInversion23
        localModeSymmetryOperations[6, 3:5, 3:5] = permutationInversion13
        localModeSymmetryOperations[6, 7:9, 7:9] = permutationInversion13
        localModeSymmetryOperations[6, 10:12, 10:12] = permutationInversion13
    end
    return localModeSymmetryOperations
end

function generateSymmetryOperationsXi(case="C3v(M)"::String)::Array{SymPy.Sym}
    permuation123::Matrix{SymPy.Sym} = SymPy.Sym[0 1 0; 0 0 1; 1 0 0]
    permuation123Dihedrals::Matrix{SymPy.Sym} = SymPy.Sym[-1/2 sqrt(3)/2; -sqrt(3)/2 -1/2]
    permuation132::Matrix{SymPy.Sym} = SymPy.Sym[0 0 1; 1 0 0; 0 1 0]
    permuation132Dihedrals::Matrix{SymPy.Sym} = SymPy.Sym[-1/2 -sqrt(3)/2; sqrt(3)/2 -1/2]
    permutationInversion12::Matrix{SymPy.Sym} = SymPy.Sym[0 1; 1 0]
    permutationInversion12Dihedrals::Matrix{SymPy.Sym} = SymPy.Sym[-1/2 sqrt(3)/2; sqrt(3)/2 1/2]
    permutationInversion23::Matrix{SymPy.Sym} = SymPy.Sym[0 1; 1 0]
    permutationInversion23Dihedrals::Matrix{SymPy.Sym} = SymPy.Sym[1 0; 0 -1]
    permutationInversion13::Matrix{SymPy.Sym} = SymPy.Sym[0 0 1; 0 1 0; 1 0 0]
    permutationInversion13Dihedrals::Matrix{SymPy.Sym} = SymPy.Sym[-1/2 -sqrt(3)/2; -sqrt(3)/2 1/2]
    if case == "C3v(M)"
        xiSymmetryOperations::Array{SymPy.Sym} = zeros(Int64, 6, 11, 11)
        for i in 1:6
            xiSymmetryOperations[i, :, :] = Matrix(1I, 11, 11)
        end
        xiSymmetryOperations[2, 3:5, 3:5] = permuation123
        xiSymmetryOperations[2, 7:9, 7:9] = permuation123
        xiSymmetryOperations[2, 10:11, 10:11] = permuation123Dihedrals
        xiSymmetryOperations[3, 3:5, 3:5] = permuation132
        xiSymmetryOperations[3, 7:9, 7:9] = permuation132
        xiSymmetryOperations[3, 10:11, 10:11] = permuation132Dihedrals
        xiSymmetryOperations[4, 3:4, 3:4] = permutationInversion12
        xiSymmetryOperations[4, 7:8, 7:8] = permutationInversion12
        xiSymmetryOperations[4, 10:11, 10:11] = permutationInversion12Dihedrals
        xiSymmetryOperations[5, 4:5, 4:5] = permutationInversion23
        xiSymmetryOperations[5, 8:9, 8:9] = permutationInversion23
        xiSymmetryOperations[5, 10:11, 10:11] = permutationInversion23Dihedrals
        xiSymmetryOperations[6, 3:5, 3:5] = permutationInversion13
        xiSymmetryOperations[6, 7:9, 7:9] = permutationInversion13
        xiSymmetryOperations[6, 10:11, 10:11] = permutationInversion13Dihedrals
    end
    return xiSymmetryOperations
end

function generateSymmetryOperationsCartesianRepresentation(zMatrix::DataFrame, case="C3v(M)"::String)
    if case == "C3v(M)"
        cartesianSymmetryOperations = zeros(Int64, 6, 18, 18)
        cartesianSymmetryOperations[1, :, :] = matrix(1I, 18, 18)
    end
    return cartesianSymmetryOperations
end

function generateXiCoordinates(localModeCoordinates::Vector{SymPy.Sym})::Vector{SymPy.Sym}
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
    
    xi::Vector{SymPy.Sym} = Vector{SymPy.Sym}()
    
    # Stretches
    push!(xi, 1.00 - exp(-a1*(localModeCoordinates[1] - rCOeq)))
    push!(xi, 1.00 - exp(-a2*(localModeCoordinates[2] - rOHeq)))
    push!(xi, 1.00 - exp(-b1*(localModeCoordinates[3] - rH1eq)))
    push!(xi, 1.00 - exp(-b1*(localModeCoordinates[4] - rH1eq)))
    push!(xi, 1.00 - exp(-b1*(localModeCoordinates[5] - rH1eq)))

    # Angles
    push!(xi, localModeCoordinates[6] - alphaCOHeq)
    push!(xi, localModeCoordinates[7] - alphaOCHeq)
    push!(xi, localModeCoordinates[8] - alphaOCHeq)
    push!(xi, localModeCoordinates[9] - alphaOCHeq)

    # Define symmeterised dihedrals
    d12::SymPy.Sym = localModeCoordinates[11] - localModeCoordinates[10]
    d23::SymPy.Sym = localModeCoordinates[12] - localModeCoordinates[11]
    d13::SymPy.Sym = localModeCoordinates[10] - localModeCoordinates[12]
    push!(xi, (2*d23 - d13 - d12)/sqrt(6))
    push!(xi, (d13 - d12)/sqrt(2))
    
    # Torsion angle
    tau::SymPy.Sym = (localModeCoordinates[10] + localModeCoordinates[11] + localModeCoordinates[12])/3
    push!(xi, tau)

    return xi
end

function generateInitialPotentialParameters(maxOrder=6::Int64, maxMultiMode=2::Int64)::DataFrame
    label::Vector{SymPy.Sym} = Vector{SymPy.Sym}()
    iPower::Vector{Int64} = Vector{Int64}()
    jPower::Vector{Int64} = Vector{Int64}()
    kPower::Vector{Int64} = Vector{Int64}()
    lPower::Vector{Int64} = Vector{Int64}()
    mPower::Vector{Int64} = Vector{Int64}()
    nPower::Vector{Int64} = Vector{Int64}()
    oPower::Vector{Int64} = Vector{Int64}()
    pPower::Vector{Int64} = Vector{Int64}()
    qPower::Vector{Int64} = Vector{Int64}()
    rPower::Vector{Int64} = Vector{Int64}()
    sPower::Vector{Int64} = Vector{Int64}()
    tPower::Vector{Int64} = Vector{Int64}()
    parameters::Vector{Float64} = Vector{Float64}()
    for i in 0:maxOrder
        for j in 0:maxOrder-i
            for k in 0:maxOrder-(i+j)
                for l in 0:maxOrder-(i+j+k)
                    for m in 0:maxOrder-(i+j+k+l)
                        for n in 0:maxOrder-(i+j+k+l+m)
                            for o in 0:maxOrder-(i+j+k+l+m+n)
                                for p in 0:maxOrder-(i+j+k+l+m+n+o)
                                    for q in 0:maxOrder-(i+j+k+l+m+n+o+p)
                                        for r in 0:maxOrder-(i+j+k+l+m+n+o+p+q)
                                            for s in 0:maxOrder-(i+j+k+l+m+n+o+p+q+r)
                                                for t in 0:maxOrder-(i+j+k+l+m+n+o+p+q+r+s)
                                                    multiMode = 0
                                                    if i >= 1
                                                        multiMode = multiMode + 1
                                                    end
                                                    if j >= 1
                                                        multiMode = multiMode + 1
                                                    end
                                                    if k >= 1
                                                        multiMode = multiMode + 1
                                                    end
                                                    if l >= 1
                                                        multiMode = multiMode + 1
                                                    end
                                                    if m >= 1
                                                        multiMode = multiMode + 1
                                                    end
                                                    if n >= 1
                                                        multiMode = multiMode + 1
                                                    end
                                                    if o >= 1
                                                        multiMode = multiMode + 1
                                                    end
                                                    if p >= 1
                                                        multiMode = multiMode + 1
                                                    end
                                                    if q >= 1
                                                        multiMode = multiMode + 1
                                                    end
                                                    if r >= 1
                                                        multiMode = multiMode + 1
                                                    end
                                                    if s >= 1
                                                        multiMode = multiMode + 1
                                                    end
                                                    if multiMode <= maxMultiMode
                                                        push!(label, SymPy.Sym("f$(i)$(j)$(k)$(l)$(m)$(n)$(o)$(p)$(q)$(r)$(s)$(t)"))
                                                        push!(iPower, i)
                                                        push!(jPower, j)
                                                        push!(kPower, k)
                                                        push!(lPower, l)
                                                        push!(mPower, m)
                                                        push!(nPower, n)
                                                        push!(oPower, o)
                                                        push!(pPower, p)
                                                        push!(qPower, q)
                                                        push!(rPower, r)
                                                        push!(sPower, s)
                                                        push!(tPower, t)
                                                        push!(parameters, 0.0000)
                                                        push!(label, SymPy.Sym("h$(i)$(j)$(k)$(l)$(m)$(n)$(o)$(p)$(q)$(r)$(s)$(t)"))
                                                        push!(iPower, i)
                                                        push!(jPower, j)
                                                        push!(kPower, k)
                                                        push!(lPower, l)
                                                        push!(mPower, m)
                                                        push!(nPower, n)
                                                        push!(oPower, o)
                                                        push!(pPower, p)
                                                        push!(qPower, q)
                                                        push!(rPower, r)
                                                        push!(sPower, s)
                                                        push!(tPower, t)
                                                        push!(parameters, 0.0000)
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    potentialParameters::DataFrame = DataFrame(Labels=label, iIndex=iPower, jIndex=jPower, kIndex=kPower, lIndex=lPower, mIndex=mPower,
                        nIndex=nPower, oIndex=oPower, pIndex=pPower, qIndex=qPower, rIndex=rPower, sIndex=sPower, tIndex=tPower, Parameters=parameters
    )
    return potentialParameters
end

function obtainTransformedPotentialTermsLocalMode(potentialParameters::DataFrame, localModeCoordinates::Vector{SymPy.Sym}, symmetryOperations::Array{Int64})::DataFrame
    totalNumberOfParameters::Int64 = size(potentialParameters)[1]
    numberOfModes::Int64 = size(localModeCoordinates)[1]
    numberOfTransformations::Int64 = size(symmetryOperations)[1] + 1 # Also counts untransformed coordinates

    # Defining xi before and after each of the symmetry operations is applied to it
    xiMatrix::Matrix{SymPy.Sym} = Matrix{SymPy.Sym}(zeros(numberOfTransformations, numberOfModes))
    xiMatrix[1, :] = generateXiCoordinates(localModeCoordinates)
    xiMatrix[2, :] = generateXiCoordinates(symmetryOperations[1, :, :]*localModeCoordinates)
    xiMatrix[3, :] = generateXiCoordinates(symmetryOperations[2, :, :]*localModeCoordinates)
    xiMatrix[4, :] = generateXiCoordinates(symmetryOperations[3, :, :]*localModeCoordinates)
    xiMatrix[5, :] = generateXiCoordinates(symmetryOperations[4, :, :]*localModeCoordinates)
    xiMatrix[6, :] = generateXiCoordinates(symmetryOperations[5, :, :]*localModeCoordinates)
    xiMatrix[7, :] = generateXiCoordinates(symmetryOperations[6, :, :]*localModeCoordinates)
    xiPowers::Matrix{SymPy.Sym} = Matrix{SymPy.Sym}(ones(totalNumberOfParameters, numberOfTransformations))
    for row in 1:totalNumberOfParameters
        if "$(string(potentialParameters[row, :Labels][1])[1])" == "f"
            for transformation in 1:numberOfTransformations
                xiPower::SymPy.Sym = 1.0
                for mode in 1:numberOfModes-1
                    xiPower = xiPower*xiMatrix[transformation, mode]^potentialParameters[row, names(potentialParameters)[1 + mode]]
                end
                xiPower = xiPower*cos(xiMatrix[transformation, 12]*potentialParameters[row, :tIndex])
                xiPowers[row, transformation] = xiPower
            end
        else
            for transformation in 1:numberOfTransformations
                xiPower::SymPy.Sym = 1.0
                for mode in 1:numberOfModes-1
                    xiPower = xiPower*xiMatrix[transformation, mode]^potentialParameters[row, names(potentialParameters)[1 + mode]]
                end
                xiPower = xiPower*sin(xiMatrix[transformation, 12]*potentialParameters[row, :tIndex])
                xiPowers[row, transformation] = xiPower            
            end
        end
    end
    potentialTerms::DataFrame = potentialParameters
    insertcols!(potentialTerms, 15, :xi => xiPowers[:, 1])
    insertcols!(potentialTerms, 16, :xi1 => xiPowers[:, 2])
    insertcols!(potentialTerms, 17, :xi2 => xiPowers[:, 3])
    insertcols!(potentialTerms, 18, :xi3 => xiPowers[:, 4])
    insertcols!(potentialTerms, 19, :xi4 => xiPowers[:, 5])
    insertcols!(potentialTerms, 20, :xi5 => xiPowers[:, 6])
    insertcols!(potentialTerms, 21, :xi6 => xiPowers[:, 7])
    return potentialTerms
end

function setupTensorFormPotential(numberOfModes=12::Int64, maxMultiMode=2::Int64, maxFourierOrder=20::Int64)::Tuple{Dict{Int64, Array{SymPy.Sym}}, Dict{Int, Vector{String}}}
    coefficients::Dict{Int64, Array{SymPy.Sym}} = Dict{Int64, Array{SymPy.Sym}}()
    firstOrderCoefficients::Array{SymPy.Sym} = zeros(2*maxFourierOrder + 1, numberOfModes-1)
    powers::Vector{Int64} = zeros(numberOfModes-1)
    label::String = "f_"
    labels::Vector{String} = Vector{String}()
    coefficientLabels::Dict{Int64, Vector{String}} = Dict{Int64, Vector{String}}() 
    for i in 1:numberOfModes-1
        powers[i] = powers[i] + 1
        for n in -maxFourierOrder:maxFourierOrder
            label = "f_"
            for j in 1:numberOfModes-1
                label = label*"$(powers[j])"
            end
            label = label*"t$(n)"
            push!(labels, label)
            firstOrderCoefficients[n+maxFourierOrder+1, i] = SymPy.Sym(label)
        end
        powers = zeros(numberOfModes)
    end
    coefficientLabels[1] = unique(labels)
    coefficients[1] = firstOrderCoefficients

    labels = Vector{String}()
    # secondOrderCoefficients::Array{SymPy.Sym} = zeros(numberOfModes, numberOfModes)
    # for i in 1:numberOfModes
    #     for j in 1:numberOfModes
    #         label = fourierType
    #         powers[i] = powers[i] + 1
    #         powers[j] = powers[j] + 1
    #         for k in 1:numberOfModes
    #             label = label*"$(powers[k])"
    #             push!(labels, label)
    #         end
    #         secondOrderCoefficients[i, j] = SymPy.Sym(label)
    #         powers = zeros(numberOfModes)
    #     end
    # end
    # coefficients[2] = secondOrderCoefficient

    return coefficients, coefficientLabels
end

function defineTorsionSymmetryOperations(torsionPower::Int64, symmetryOperations::Array{SymPy.Sym})::Array{SymPy.Sym}
    fullSymmetryOperations::Array{SymPy.Sym} = symmetryOperations
    fullSymmetryOperations[2, :, :] = fullSymmetryOperations[2, :, :]*exp(2*pi*torsionPower*1im/3)
    fullSymmetryOperations[3, :, :] = fullSymmetryOperations[3, :, :]*exp(-2*pi*torsionPower*1im/3)
    return fullSymmetryOperations
end

function solveForSymmetricCoefficients(expansionCoefficientsTensorForm::Dict{Int64, Array{SymPy.Sym}}, expansionCoefficients::Dict{Int64, Vector{String}}, symmetryOperations::Array{SymPy.Sym}, numberOfModes=12::Int64)
    maxFourierOrder::Int64 = (size(expansionCoefficientsTensorForm[1])[1] - 1) / 2
    numberOfSymmetryOperations::Int64 = size(symmetryOperations)[1]
    equationsToSolve::Dict{Int64, Array{SymPy.Sym}} = Dict{Int64, Array{SymPy.Sym}}()
    firstOrderEquations::Array{SymPy.Sym} = zeros(numberOfSymmetryOperations, 2*maxFourierOrder + 1, numberOfModes-1)
    for n in -maxFourierOrder:maxFourierOrder
        # fullSymmetryOperations::Array{SymPy.Sym} = defineTorsionSymmetryOperations(n, symmetryOperations)
        for operationNumber in 1:numberOfSymmetryOperations
            firstOrderEquations[operationNumber, n+maxFourierOrder+1, :] = reshape(expansionCoefficientsTensorForm[1][n+maxFourierOrder+1, :], 1 ,11) - reshape(expansionCoefficientsTensorForm[1][n+maxFourierOrder+1, :], 1, 11)*symmetryOperations[operationNumber, :, :]
        end
    end
    println(symmetryOperations[1, :, :])
    println(symmetryOperations[2, :, :])
    println(symmetryOperations[2])
    println("Der fugl ist sehr klug")
    for equation in firstOrderEquations
        println(typeof(equation))
        println(equation)
    end
    firstOrderSolutions = solve(firstOrderEquations, SymPy.Sym.(expansionCoefficients[1]))
    print("Die katze ist sehr gutt")
    println(firstOrderSolutions)
end