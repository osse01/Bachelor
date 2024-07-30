using Plots
include("laverySplinesBi.jl")

len = 4;
xdata = [0,1,2,3]
ydata = [0,1,2,3]

zdata = [0 0 0 0;
         0 1 1 0;
         0 1 1 0;
         0 0 0 0]
biCubicSpline(xdata, ydata, zdata, 20, 0.01)
# Create the plot
#gr()
#p = plot(xdata, ydata, zdata, seriestype=:scatter, label="Data points")
#plot!(xs, ys, label="Lavery spline")
#xaxis!("x")  # Set x-axis label
#yaxis!("y")  # Set y-axis label
#zaxis!("z")  # Set z-axis label
#title!("Points interpolated with Bi-variate Lavery Splines")
#display(p)
#savefig("3D_data.png")