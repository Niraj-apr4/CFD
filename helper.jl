 function init_temps!(nodes)
    """
    Initialize internal and Boundary 4 nodes temperature 
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

function calculate_internal_temp(P_node_loc)
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

# equations for Boundary 4 
function calculate_B4_temp(P_node_loc)
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

function calculate_temp!(nodes)
    # iterate only on B4(Newmann) and internal nodes
    for i = 1:n-2 , j = 2:n-2 
        current_node = nodes[i,j] 

        if i == 1 # check for Boundary 4 
            new_temp = calculate_B4_temp([i,j])

        else # internal nodes 
            new_temp = calculate_internal_temp([i,j])
        end

        # update the current node
        current_node[3] = new_temp 
    end
end


function experiment!(tolerance)
    # needed global variables
    global n

    # n random nodes are selected from [1, n-1^2] nodes
    test_nodes_id = rand(1:(n-1)^2 , n)

    # initialize and undef array for comparison
    init_nodes_temp = Array{Float64}(undef,n)  
    updated_nodes_temp = Array{Float64}(undef,n)  

    # setup the loop counter for plotting
    loop_counter = 0 

    while true  
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
