using Plots
include("laverySplinesBi.jl")
include("testData.jl")

xdata = xdata_2
ydata = ydata_2
zdata = zdata_2

N = 20
M = N
spline = biCubicSpline(xdata, ydata, zdata, N, 0.0001)

x = range(xdata[1], xdata[end], N)
y = range(ydata[1], ydata[end], M)
z = evaluate(spline, N, M)

# Create the plot
#gr()
xdata = xPlotData_2
ydata = yPlotData_2
zdata = zPlotData_2

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