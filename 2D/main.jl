using Plots
include("laverySplinesBi.jl")

len = 400;
xdata = range(0, stop=15, length=len);
ydata = range(0, stop=15, length=len);

zdata = cos.(xdata) .+ sin.(ydata);

# Create the plot
#gr()
p = plot(xdata, ydata, zdata, seriestype=:scatter, label="Data points")
#plot!(xs, ys, label="Lavery spline")
xaxis!("x")  # Set x-axis label
yaxis!("y")  # Set y-axis label
zaxis!("z")  # Set z-axis label
title!("Points interpolated with Bi-variate Lavery Splines")
display(p)
savefig("3D_data.png")