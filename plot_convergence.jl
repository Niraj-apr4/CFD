"""
cfd2.jl
objectives : for calculating the number of iteration
             required for given tolerance and plot in logarithmic
             scale
"""
include("cfd2.jl")

tolerance_array = [0.1 0.01 0.001 0.0001 0.00001 0.000001 0.0000001 0.00000001 1e-9 1e-10 1e-11]
loops_array = similar(tolerance_array)

for i in 1: length(tolerance_array)
    loops_array[i] = experiment(tolerance_array[i])
end

function plot_convergence(loops_array, tolerance_array)
    plot(loops_array, tolerance_array,
         xlabel="Iteration Number",
         ylabel="Convergence Tolerance",
         title="Convergence History(logarithmic plot)",
         yscale=:log10,                    # Log scale for tolerance
         color=:blue,
         marker=:square,
         grid=:true,
         framestyle=:box,
         markersize=5,
         size=(800, 600),
         fontfamily="Computer Modern",
         titlefontsize=16,
         guidefontsize=14,
         tickfontsize=12,
         legend=false)
    
    # annotation
    hline!([tolerance_array[7]], linestyle=:dash, color=:green, label="Target")
    vline!([loops_array[7]], linestyle=:dash, color=:green, label="Target")
end

# Usage:
p = plot_convergence(loops_array, tolerance_array)
savefig(p, "converge_log_test.png")
display(p)
