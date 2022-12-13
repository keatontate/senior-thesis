### define lattice operators, and store glide plane translations as a fourth column since space group 39 is nonsymmorphic

C2x = [1 0 0 0
       0 -1 0 0
       0 0  -1 0]
C2z = [-1 0 0 0
       0 -1 0 0
       0 0  1 0]
C2y = [-1 0 0 0
       0 1 0 0
       0 0  -1 0]
σx = [-1 0 0 0
       0 1 0 0.5
       0 0 1 0]
σy = [1 0 0 0
      0 -1 0 0.5
      0 0 1 0]
σz = [1 0 0 0
      0 1 0 0.5
      0 0 -1 0 ]
I = [-1 0 0
      0 -1 0
      0 0 -1]
E = [1 0 0 0
     0 1 0 0
     0 0 1 0]
# conventional unit cell basis vectors
c_lVs = [8.512 0 0 
        0 9.242 0
        0 0 8.313]

# primitive unit cell basis vectors (for the fractional atom positions)
p_lVs = [0 4.621 4.621 
        0 -4.1565 4.1565
        8.3123 0 0]

# p_lVs = [8.512 0 0 
#         0 4.621 4.621
#         0 -4.1565 4.1565]

# the positions of the atoms in fractional coords with primitive lattice vectors
vecs = [0.3832 0.1366 0.1804
        0.6398 0.8862 0.2295
        0.4023 0.6165 0.3411
        0.5352 0 0.5
        0.6958 0.75 0.0613
        0.4633 0.75 0.0266
        0.4351 0.9355 0.1182
        0.4706 0.75 0.6311
        0.2759 0.75 0.1870
        0.1129 0.8972 0.1287
        0.1334 0.1075 0.6443
        0.7434 0 0
        0.3008 0.0604 0.3845
        0.2759 0.75 0.5035
        0.5782 0.75 0.3917]

# store all our symmetry operators to loop over later
syms = [C2x, σy, σz, E]

keepvecs = zeros(44,3)

i = 0
for pos in eachrow(vecs)
    println("new pos")
    j = 0
    for sym in syms
        cart = c_lVs * pos  # Get vector to cartesian coords
        newPos = sym[1:3,1:3] * cart  # Apply symmetry operator
        direct = inv(c_lVs) * newPos  # Go back to direct coords
        inCell = mod.(direct + sym[:,4],1) # Map back into cell
        #println("symmetry operator")
        #display(sym)
        #println("before rotation (fractional)")
        #display(pos)
        #println("before rotation (cartesian)")
        #display(cart)
        #println("after rotation (cartesian)")
        #display(newPos)
        #println("After rotation(fractional)")
        #display(inCell)
        
        ### use this one if you want the conventional unit cell, it includes the glide plane translate
        # inCellTwo = mod.(direct + sym[:,4] + [0.5, 0.5 , 0],1)  # Get the second translate ### Move this outside the loop
        
        ### this one is for generating the primitive cell positions
        inCellTwo = mod.(direct + sym[:,4],1)  # Get the second translate ### Move this outside the loop
        keep = !any([all(isapprox.(inCell,j)) for j in eachrow(keepvecs[1:i,:])])
        
        #println("should we keep it?")
        #println(keep)
        if keep
            j += 1
            global i += 1
            keepvecs[i,:] = inCell
        else
            println("throw it out")
        end
        
        # keepTwo = !any([all(isapprox.(inCellTwo,j)) for j in eachrow(keepvecs[1:i,:])])
        # if keepTwo
        #     i += 1
        #     keepvecs[i,:] = inCellTwo
        # else
        #     println("throw it out")
        # end
        # #keepvecs[i,:] = inCellTwo  # Add vector to list
        # #i+=1

    end
    println("kept this many:")
    println(j)
#    println("check")
#    println(i)
    #[]
    #println(size(unique(keepvecs, dims = 1)))
end

println(i)
# @assert i==88
@assert i==44

# keepvecs_cart = keepvecs * c_lVs

# # create a supercell here to compare distances

# using Distances
# R = pairwise(Euclidean(), keepvecs_cart, dims=1)

# R = R[1:i,1:i] # trim R so I don't get more than the atoms we keep

# testR = R[findall(x -> x < 1.77, R)] # this should be empty except for the trace

# ### need to implement some way to find the smallest distance off the trace of this matrix
# smallest_dist, sdist_idx = R[argmin(R)], argmin(R)
# #argmin(R), "\nMin distance: ", R[argmin(R)])

# using Plots
# #plotlyjs()
# heatmap(R, title="Shortest Distance: $smallest_dist A\n $sdist_idx", aspect_ratio=1)

### perform a change of basis to get the lattice in terms of primitive cell vectors
### Ax = By (conventional basis) (conventional positions) = (primitive basis) (primitive positions)
### y = B_inv(Ax)
pvecs = inv(p_lVs) * c_lVs * keepvecs'

# now we translate back into the unit cell again by clever mod like earlier
pvecs = mod.(pvecs,1)

### write the output to file for use in creating POSCAR
using DelimitedFiles
cd("/home/keaton/Documents/School/SPRING2022/SeniorResearch/projects/lithium-project/structures/oc88/")
writedlm("oC88.txt",pvecs')
display(pvecs[begin:25,:])
display(size(unique(pvecs,dims = 1)))
# mod(-0.1366,-1)


### these are for testing, showing how unique works and how sensitive isapprox() is

# A = [1 2 1
#      1 2 1
#      3 2 1]
# B = unique(A,dims = 1)
# display(B)

# b = [0.5352, 0.0, 0.5]
# c = [0.535200001,0.0,0.5]
# all(isapprox.(b,c))

# testx = [0.5; 
#          0.5; 
#          0.5]

# display(c_lVs * testx)
# display(p_lVs * testx)

