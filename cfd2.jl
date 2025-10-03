
"""
Objectives : COMPUTATIONS

           Assumes points and nodes are available from mesh_geometery.jl

        STEP 1 : identify the boundary nodes and apply the
                 boundary conditions  

                 Boundary 1 : T_1 = 15
                 Boundary 2 : T_2 = 10 
                 Boundary 3 : T_3 = 5(1-y/H) + 15 * sin(pi*y/H)

        STEP 2 : write the equation for boundary 4

        STEP 3 : write the equation for internal nodes

        STEP 4 : setup the conditons for tolerance
                 Approach :
                 1. pick n random nodes from grid
                 2. save temperature before each iteration
                 3. find the temperature after iteratrion 
                 4. diff = after_iteration_temperature - before_temperature
                 5. elementwise square each difference
                    diff.^2
                 6. max(diff.^2) < tolerance   
                 and

                 prepare the required helper function for computation

        STEP 5 : perform the computations 
"""

include("mesh_geometry.jl")
n = 6
L = 6 
H = 6
mesh2 = generate_2Dmesh(L,H,n)
nodes = mesh2[2]
points = mesh2[1]


# identify the boundary nodes and apply Boundary condition
Boundary1 = nodes[:   ,   1]
Boundary2 = nodes[n-1 ,   :]
Boundary3 = nodes[:   , n-1]
Boundary4 = nodes[1   ,   :]

# STEP 1 ##### process boundary values (Dirichlet Condition)

# apply the Dirichelt boundary conditions 
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

# calculation for Newmann Boundary condtions
# boundary 4

# equations for Boundary 4 
function B4Temp(P_node_loc)
    """
    """
    global n,H

    # node selection
    i = P_node_loc[1]
    j = P_node_loc[2]
    P_node = nodes[i,j] # current node

    # required current node parameter 
    y = P_node[2]
    k  = 16(y/H + 1)

    # east node parameter
    E_node = nodes[i+1,j]
    TempE = E_node[3]
    delta_xPE = E_node[1] - P_node[1]  

    if j != 1 && j != n-1 
        # required Adjacent node 
        S_node = nodes[i   , j-1]
        N_node = nodes[i   , j+1]
        # no west node since at boundary 4

        # required node parameters
        delta_ySP = P_node[2] - S_node[2]
        delta_yPN = N_node[2] - P_node[2]

        delta_y = delta_ySP/2 + delta_yPN/2
        Ae = delta_y

    elseif j == 1 
        # no south node
        N_node = nodes[i   , j+1]
        delta_yPN = N_node[2] - P_node[2]

        Ae = delta_yPN/2

    elseif j == n-1 
        # no north node
        S_node = nodes[i   , j-1]
        delta_ySP = P_node[2] - S_node[2]

        Ae = delta_ySP/2
    end

    # Flux at Boundary4 is q = -5000
    TempP = TempE - 5000 * delta_xPE/(k * Ae) 

end

# equations for internal nodes >>>
function internalTemp(P_node_loc)
    """
    # Arguments
    - `P_node_loc`: 2 element array [i,j] representing location of p_node
    # outputs
    - TempP : updated temperature at node P
    """
    # use global parameters 
    global H, nodes 

    i = P_node_loc[1] # current x cordinate
    j = P_node_loc[2] # current y cordinate 

    P_node = nodes[i,j] # current node
    W_node = nodes[i-1 , j  ] # Adjacent nodes
    S_node = nodes[i   , j-1]
    E_node = nodes[i+1 , j  ]
    N_node = nodes[i   , j+1]

    # node geometry
    delta_xWP = P_node[1] - W_node[1]
    delta_xPE = E_node[1] - P_node[1]  
    delta_ySP = P_node[2] - S_node[2]
    delta_yPN = N_node[2] - P_node[2]

    # TODO Please confirm if this is correct 
    delta_x = delta_xWP/2 + delta_xPE/2
    delta_y = delta_ySP/2 + delta_yPN/2

    Aw = delta_y # face areas in 2D case
    Ae = delta_y
    An = delta_x
    As = delta_x

    # calculate k value according to y-cordinate of nodes   
    kw  = 16*(W_node[2]/H + 1)
    ke  = 16*(E_node[2]/H + 1)
    ks  = 16*(S_node[2]/H + 1)
    kn  = 16*(N_node[2]/H + 1)

    # TODO please confirm if this is correct
    Su = 0
    # P_node[3] is current temperature at P 
    # caution devide by 0 error if previous temp is 0
    Sp = (1.5 * delta_x *  delta_y)/P_node[3]

    aW = kw * Aw/delta_xWP 
    aE = ke * Ae/delta_xPE 
    aS = ks * As/delta_ySP 
    aN = kn * An/delta_yPN 

    aP = aW + aE + aS + aN - Sp

    TempW = W_node[3] # temperature at adjacent nodes 
    TempS = S_node[3]
    TempE = E_node[3]
    TempN = N_node[3]


    # update and return  TempP 
    TempP = (aW* TempW + aE * TempE + aS * TempS +
        aN * TempN + Su) / aP
end
# <<<

# TODO
# after applying Boundary Condtion => initial Processing
# assign initial temperature at all nodes except nodes at dirichlet
# Boundary conditons

 function initialize_internal_temps!(nodes)
    """
    Initialize internal and Boundary 4 nodes with reasonable guesses
    based on linear interpolation from Dirichlet boundaries
    """
    global n, H
    
    n_rows, n_cols = size(nodes)
    
    # Get boundary values for interpolation
    T_bottom = 15.0  # Boundary 1 (constant)
    T_top = 10.0     # Boundary 3 (constant)
    
    for i = 1:n_rows, j = 1:n_cols
        # Skip already-set Dirichlet boundaries
        if i == n_rows  # Boundary 2 (already set)
            continue
        elseif j == 1   # Boundary 1 (already set)
            continue
        elseif j == n_cols  # Boundary 3 (already set)
            continue
        end
        
        # For Boundary 4 and internal nodes
        y = nodes[i,j][2]
        
        # Linear interpolation between top and bottom
        nodes[i,j][3] = T_bottom + (T_top - T_bottom) * (y / H)
    end
end   


# # STEP 4 ##### test for convergence  and required helper functions

function calculate_temp!(nodes)
    """
    Calculates Temperature at
    1. Boundary 4
    2. All internal nodes
    mutates the nodes array 
    """
    for i = 1:n-2 , j = 2:n-2 # iterated on boundary 4 and internal nodes 
        current_node = nodes[i,j] 

        if i == 1 # check for Boundary 4 
            new_temp = B4Temp([i,j])

        else # internal nodes 
            new_temp = internalTemp([i,j])
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
    # Approach
     Approach is to set a tolerance value and
    continue the iteration untill the temperature at nodes
    start converging ,

    - After each iteration test at random nodes if the difference
    of convergence is with in an acceptable range defined by
    parameter tolerance  

    - uses calculate_temp! which mutates the original nodes 

    # Arguments
    -`tolerance::Float`: set up the tolerance

    # Returns and Operations

    - call calculate_temp! function which mutates the nodes array
    with new temperatures and continues till the value start to converge
    at randomly selected n nodes

    -`loop_counter::Int` : return the number of loops required

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
tolerance = 0.00000001
# start the computations
experiment!(tolerance)
