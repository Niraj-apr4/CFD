include("mesh_geometry.jl")
using Plots

mesh1 = generate_2Dmesh(6,6,6)

points1 = mesh1[1]
nodes1 = mesh1[2]
