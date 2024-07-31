using Plots
include("laverySplinesBi.jl")
include("testData.jl")

xdata = xdata_1
ydata = ydata_1
zdata = zdata_1

# len = 4;
# xdata = [0,1,2,3,4,5,6]
# ydata = [0,1,2,3,4,5,6]
# zdata = [0 0 0 0 0 0 0;
#          0 0 0 0 0 0 0;
#          0 0 1 1 1 0 0;
#          0 0 1 1 1 0 0;
#          0 0 1 1 1 0 0;
#          0 0 0 0 0 0 0;
#          0 0 0 0 0 0 0]
N = 20
M = N
spline = biCubicSpline(xdata, ydata, zdata, N, 0.0001)
bx = spline.bx
by = spline.by
x = range(xdata[1], xdata[end], N)
y = range(ydata[1], ydata[end], M)
z = evaluate(spline, N, M)

# Create the plot
#gr()
# xdata = [0,1,2,3,4,5,6, 0,1,2,3,4,5,6, 0,1,2,3,4,5,6, 0,1,2,3,4,5,6, 0,1,2,3,4,5,6, 0,1,2,3,4,5,6, 0,1,2,3,4,5,6]
# ydata = [0,0,0,0,0,0,0, 1,1,1,1,1,1,1, 2,2,2,2,2,2,2, 3,3,3,3,3,3,3, 4,4,4,4,4,4,4, 5,5,5,5,5,5,5, 6,6,6,6,6,6,6]
# zdata = vec(zdata)
xdata = xPlotData_1
ydata = yPlotData_1
zdata = zPlotData_1

# , quiver = (vec(bx), vec(by), ones(1,len))
plot(x, y, z, seriestype=:surface, label="Lavery spline")
plot!(xdata, ydata, zdata, seriestype=:scatter, label="data points", camera=(0,90))
xaxis!("x")  # Set x-axis label
yaxis!("y")  # Set y-axis label
zaxis!("z")  # Set z-axis label
title!("Points interpolated with Bi-variate Lavery Splines")
#display(p)
savefig("bi_lavery_top.png")

plot(x, y, z, seriestype=:surface, label="Lavery spline")
plot!(xdata, ydata, zdata, seriestype=:scatter, label="data points")
xaxis!("x")  # Set x-axis label
yaxis!("y")  # Set y-axis label
zaxis!("z")  # Set z-axis label
title!("Points interpolated with Bi-variate Lavery Splines")
#display(p)
savefig("bi_lavery_angle.png")