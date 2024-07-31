using Plots
include("laverySplinesBi.jl")

len = 4;
xdata = [0,1,2,3]
ydata = [0,1,2,3]

zdata = [0 0 0 0;
         0 1 1 0;
         0 1 1 0;
         0 0 0 0]
N = 20
M = N
spline = biCubicSpline(xdata, ydata, zdata, N, 0.01)
x = range(xdata[1], xdata[end], N)
y = range(ydata[1], ydata[end], M)
z = evaluate(spline, N, M)

# Create the plot
#gr()
xdata = [0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]
ydata = [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3]
zdata = vec(zdata)
plot(x, y, z, seriestype=:surface, label="Lavery spline")
p = plot!(xdata, ydata, zdata, seriestype=:scatter, label="Data points")
xaxis!("x")  # Set x-axis label
yaxis!("y")  # Set y-axis label
zaxis!("z")  # Set z-axis label
title!("Points interpolated with Bi-variate Lavery Splines")
display(p)
savefig("bi_lavery.png")