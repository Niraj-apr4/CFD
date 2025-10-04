
include("contour_plot.jl")
include("mesh_geometry.jl")
include("helper.jl")

######################
# STEP 1 create mesh # 
######################
n = 30 
L = 1 
H = 0.5
sx = 1
sy = 2 

mesh2 = generate_2Dmesh(L,H,n,sx,sy)
nodes = mesh2[2]
points = mesh2[1]

# plot_mesh(nodes)
# contour_plot(nodes)

#####################


##########################################
# STEP 2 identify the boundary nodes and #
# apply  Dirichelt Boundary condition    #
##########################################
# identify the boundary 
Boundary1 = nodes[:   ,   1] 
Boundary2 = nodes[n-1 ,   :]
Boundary3 = nodes[:   , n-1]
Boundary4 = nodes[1   ,   :]

# Apply Dirichlet Condition
# Boundary 1
for i = 1:length(Boundary1)
    #set Boundary1 temp value to 15
    Boundary1[i][3] = 15
end

# boundary 2
for i = 1:length(Boundary2)
    global H
    y = nodes[i][2] 
    # Boundary2 temp is function of y 
    Boundary2[i][3] = 5*(1-y/H) + 15 * sin(pi*y/H)
end

# boundary 3
for i = 1:length(Boundary3)
    #fix Boundary3 temp value to 10
    Boundary3[i][3] = 10
end

# after applying dirichlet conditions
# contour_plot(nodes)

############################


########################################################
# STEP 3 initialize temperature at all remaining nodes #
########################################################

# Approach : linear approximation is applied
# from Boundary1 to Boundary3
# function : init_temps!(nodes) (file: helper.jl)

init_temps!(nodes)
# demonstrate the initial temperature distributions
# contour_plot(nodes)

##############################################


#############################################################
# STEP 4 Perform calculation until it  converges with given #
# tolerance                                                 #
#############################################################

# set the tolerance here 
tolerance = 0.00001

# start the computations
experiment!(tolerance)

# contour_plot(nodes)



##########################################################
# STEP 5 setup different tolerance and measure the no of #
# loops it takes to converge                             # 
##########################################################

# tolerance_array = [0.1 0.01 0.001 0.0001 0.00001 0.000001 0.0000001 0.00000001 1e-9 1e-10 1e-11]

# TODO something strange is happening here fix
# tolerance_array = [0.1 0.01 0.001 0.0001 0.00001 0.000000001]
# loops_array = similar(tolerance_array)

# for i in 1: length(tolerance_array)
#     loops_array[i] = experiment!(tolerance_array[i])
# end
