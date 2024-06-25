function computeCartesianCoordinates(localModeCoordinates, case="bond-fixed-dihedral")
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

function generateSymmetryOperationsLocalModeRepresentation(case="C3v(M)")
    permuation123 = Int64[0 1 0; 0 0 1; 1 0 0]
    permuation132 = Int64[0 0 1; 1 0 0; 0 1 0]
    permutationInversion12 = Int64[0 1; 1 0]
    permutationInversion23 = Int64[0 1; 1 0]
    permutationInversion13 = Int64[0 0 1; 0 1 0; 1 0 0]
    if case == "C3v(M)"
        localModeSymmetryOperations = zeros(Int64, 6, 12, 12)
        localModeSymmetryOperations[1, :, :] = Matrix(1I, 12, 12)
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

function generateSymmetryOperationsCartesianRepresentation(zMatrix, case="C3v(M)")
    if case == "C3v(M)"
        cartesianSymmetryOperations = zeros(Int64, 6, 18, 18)
        cartesianSymmetryOperations[1, :, :] = matrix(1I, 18, 18)
    end
    return cartesianSymmetryOperations
end

function generateInitialPotentialParameters(maxOrder=6, maxMultiMode=2)
    label = String[]
    iPower = Vector{Int64}()
    jPower = Vector{Int64}()
    kPower = Vector{Int64}()
    lPower = Vector{Int64}()
    mPower = Vector{Int64}()
    nPower = Vector{Int64}()
    oPower = Vector{Int64}()
    pPower = Vector{Int64}()
    qPower = Vector{Int64}()
    rPower = Vector{Int64}()
    sPower = Vector{Int64}()
    tPower = Vector{Int64}()
    parameters = Vector{Float64}()
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
                                                    if t >= 1
                                                        multiMode = multiMode + 1
                                                    end
                                                    if multiMode <= maxMultiMode
                                                        push!(label, "f$(i)$(j)$(k)$(l)$(m)$(n)$(o)$(p)$(q)$(r)$(s)$(t)")
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
                                                        push!(label, "h$(i)$(j)$(k)$(l)$(m)$(n)$(o)$(p)$(q)$(r)$(s)$(t)")
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
    potentialParameters = DataFrame(Labels=label, i=iPower, j=jPower, k=kPower, l=lPower, m=mPower,
                        n=nPower, o=oPower, p=pPower, q=qPower, r=rPower, s=sPower, t=tPower, Parameters=parameters
    )
    return potentialParameters
end