# stretchesSpacing::Vector{Float64} = [0.0,-0.05,0.05,-0.1,0.10,-0.20,0.20,-0.30,0.30,0.40,0.50,0.70]
# angleSpacing::Vector{Float64} = [0.0,5.00,-5.0,10.0,-10.0,20.0,-20.0,30.0,-30.0,50.0,-50.0,60.0]
# dihedralSpacing::Vector{Float64} = [0.0,5.0,-5.0,10.0,-10.0,40.0,-40.0,60.0,-60.0,80.0,-80.0,100.0]
# torsionSpacing::Vector{Float64} = [0.00, 10.00,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.,110.0,120.0]

stretchesSpacing::Vector{String} = ["0.0000", "0.0100", "-0.0100", "-0.0200", "0.0200", "0.0300", 
    "-0.0300", "0.0400", "-0.0400", "0.0500", "-0.0500", "0.0600", "0.0700"]
angleSpacing::Vector{String} = ["0.0000", "1.0000", "-1.0000", "2.0000", "-2.0000", "3.0000", 
    "-3.0000", "4.0000", "-4.0000", "5.0000", "-5.0000", "6.0000", "-6.0000", "7.0000", 
    "-7.0000", "8.0000", "-8.0000", "9.0000", "-9.0000", "10.0000", "-10.0000"]
dihedralSpacing::Vector{String} = ["0.0000", "1.0000", "-1.0000", "2.0000", "-2.0000", "3.0000", 
    "-3.0000", "4.0000", "-4.0000", "5.0000", "-5.0000", "6.0000", "-6.0000", "7.0000", 
    "-7.0000", "8.0000", "-8.0000", "9.0000", "-9.0000", "10.0000", "-10.0000"]

stretchesGrid::Int64 = size(stretchesSpacing)[1]
angleGrid::Int64 = size(angleSpacing)[1]
for i in 1:stretchesGrid
    for j in 1:stretchesGrid
        for k in 1:stretchesGrid
            for l in 1:stretchesGrid
                for m in 1:angleGrid
                    for n in 1:angleGrid
                        for o in 1:angleGrid
                            for p in 1:angleGrid
                                for q in 1:angleGrid
                                    println("$(stretchesSpacing[i])  $(stretchesSpacing[j])  $(stretchesSpacing[k])  $(stretchesSpacing[l])  $(stretchesSpacing[k])  $(angleSpacing[m])  $(angleSpacing[n])  $(angleSpacing[o])  $(angleSpacing[n])  $(angleSpacing[p])  $(angleSpacing[q])")
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end