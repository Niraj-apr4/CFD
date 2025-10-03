"""
contour_plot.jl
objectives : for plotting contour plot with different
            mesh geometry and precision
"""

include("cfd2.jl")

# IMPORTANT! please check step 5 of cfd2.jl file 
# before proceeding onwards

# delta_x,delta_y ,n and tolerance are defined
# in mesh_geometry.jl file  please redefine there
# to get the desired result

# redefined again 
tolerance = 0.00001

# run experiment 
experiment!(tolerance)
    
function contourPlot(data)
    # Get unique coordinates
    x_vals = sort(unique([p[1] for p in data]))
    y_vals = sort(unique([p[2] for p in data]))
    
    # Create temperature matrix
    temp_matrix = zeros(length(y_vals), length(x_vals))
    for point in data
        i = findfirst(y -> y ≈ point[2], y_vals)
        j = findfirst(x -> x ≈ point[1], x_vals)
        temp_matrix[i, j] = point[3]
    end
    
    # Professional contour plot
    c = contourf(x_vals, y_vals, temp_matrix,
            title="Temperature Distribution",
            xlabel="X Coordinate (m)", 
            ylabel="Y Coordinate (m)",
            colorbar_title="Temperature (°C)",
            aspect_ratio=:equal,
            levels=25,
            color=:hot,
            size=(700, 600),
            dpi=300,
            linewidth=0,
            fontfamily="Computer Modern")
    savefig(c,"plot1.png")
end

# Usage:
contourPlot(nodes)
