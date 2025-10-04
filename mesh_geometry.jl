
"""
mesh_geometry.jl  

Objective : Define  the  Mesh Geomety
    STEP 1 : first create a differential 2D Control Volume with 

    n = required number of such differential control volumes
    required to construct full CV

    STEP 2 : Compute the computational nodes for each
             differential control volume
             computational nodes are assigned at the mid point
"""

function generate_2Dmesh(L,H,n,sx=1,sy=1)
    """
    # Arguments
    - `L::Float`: length of full cv
    - `H::Float`: Height of full cv
    - `n::Integer`: number of differential cv
                    more n more  finer the mesh is 
    # Outputs
    [points,nodes]
    - `points::Array{Array{Float64}}`:
    - `nodes::Array{Array{Float64}}`:

    """

    # STEP 1 ##### 

    # initialize An undefined array called "points"
    # with size  nxn 
    # and its element is itself an array of type {Array{Float64}}
    # with following characteristics
    #      element[1] = x_cordinate
    #      element[2] = y_cordinate
    #      element[3] = Temperature (initially set to 0)
    # 
    points = Array{Array{Float64}}(undef,n,n)

    # create the CV grid >>> 
    if sx == 1
        #unifrom
        x_array = range(0, L, n)  
    else
        # Geometric stretching
        x_array = [L * (sx^(i-1) - 1) / (sx^(n-1) - 1) for i in 1:n]
    end
    
    if sy == 1
        #uniform
        y_array = range(0, H, n) 
    else
        # Geometric stretching
        y_array = [H * (sy^(j-1) - 1) / (sy^(n-1) - 1) for j in 1:n]
    end

    x_array = range(0,L,n)
    y_array = range(0,H,n)
    for i = 1:n 
        for j = 1:n 
            points[i,j] = [x_array[i]*sx, y_array[j]*sy ,0]
        end
    end
    # now the points array consists of CV grids with initial temp
    # assigned to be  0 

    # STEP 2 ##### 

    # Compute  computational nodes for each cv >>>
    nodes = Array{Array{Float64}}(undef,n-1,n-1)
    for i = 1:n-1 ,j = 1:n-1
        init_x = points[i,j][1]
        init_y = points[i,j][2]

        final_x = points[i+1,j][1]
        final_y = points[i,j+1][2]

        delta_x = final_x - init_x
        delta_y = final_y - init_y

        nodes[i,j] = [points[i,j][1] + delta_x/2 , points[i,j][2] +
            delta_y/2 , 0]
    end
    # <<<
    return[points,nodes]
end
