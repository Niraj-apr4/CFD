
include("contour_plot.jl")
include("mesh_geometry.jl")
include("helper.jl")

# STEP 1 create mesh
n = 6 
L = 6 
H = 6
mesh2 = generate_2Dmesh(L,H,n)
nodes = mesh2[2]
points = mesh2[1]

# plot_mesh(nodes)

# contour_plot(nodes)


# STEP 2 identify the boundary nodes and
# apply  Dirichelt Boundary condition
Boundary1 = nodes[:   ,   1]
Boundary2 = nodes[n-1 ,   :]
Boundary3 = nodes[:   , n-1]
Boundary4 = nodes[1   ,   :]

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


# boundary 4 (Newmann Condition)
# Newmann condition which is defined in helper.jl
# function is written to iterate with given flux 
# function name : calculate_B4_Temp (file: helper.jl)


# STEP 3 initialize temperature at all remaining nodes

# Approach : linear approximation is applied
# from Boundary1 to Boundary3
# function : init_temps!(nodes) (file: helper.jl)

init_temps!(nodes)
# contour_plot(nodes)

# STEP 4 Perform calculation until it  converges with given
# tolerance 

function calculate_temp!(nodes)
    """
    Calculates Temperature at
    1. Boundary 4 (Newmann condition)
    2. All internal nodes
    mutates the nodes array 
    """
    # iterate only on B4(Newmann) and internal nodes
    for i = 1:n-2 , j = 2:n-2 
        current_node = nodes[i,j] 

        if i == 1 # check for Boundary 4 
            new_temp = calculate_B4_temp([i,j])
            # function :calculate_B4_temp (file: helper.jl)

        else # internal nodes 
            new_temp = calculate_internal_temp([i,j])
            # function :calculate_internal_temp (file: helper.jl)
        end

        # update the current node
        current_node[3] = new_temp 
    end
end

# n random nodes are selected from [1, n-1^2] nodes
# for purpose of testing convergence 
test_nodes_id = rand(1:(n-1)^2 , n)
init_nodes_temp = Array{Float64}(undef,n)  
updated_nodes_temp = Array{Float64}(undef,n)  

function experiment!(tolerance)
    """
    """
    # setup the loop counter for plotting
    loop_counter = 0 

    while true  
        global init_nodes_temp, updated_nodes_temp
        ## before each iterations 
        for i = 1:length(test_nodes_id) 
            init_nodes_temp[i] = nodes[test_nodes_id[i]][3]
        end

        ###################
        calculate_temp!(nodes)
        ###################

        # after iteration temp is updated
        for i = 1:length(test_nodes_id) 
            updated_nodes_temp[i] = nodes[test_nodes_id[i]][3]
        end

        test = abs(maximum((updated_nodes_temp -
            init_nodes_temp).^2)) < tolerance

        @show test

        if  test
            break
        end

        init_nodes_temp = updated_nodes_temp[:]
        loop_counter += 1
    end
    @show loop_counter
    loop_counter
end

# STEP 5 ##### Computations

# set the tolerance here 
tolerance = 0.0001
# start the computations
experiment!(tolerance)

# contour_plot(nodes)

# tolerance_array = [0.1 0.01 0.001 0.0001 0.00001 0.000001 0.0000001 0.00000001 1e-9 1e-10 1e-11]

# tolerance_array = [0.1 0.01 0.001 0.0001 0.00001]
# loops_array = similar(tolerance_array)

# for i in 1: length(tolerance_array)
#     loops_array[i] = experiment!(tolerance_array[i])
# end
