using Plots
include("laverySplinesBi.jl")
include("testData.jl")

xdata = xdata_1
ydata = ydata_1
zdata = zdata_1
lenX = length(xdata)
lenY = length(ydata)

N = 100
M = N
#spline = Data(xdata, ydata, zdata, zeros(lenX,lenY), zeros(lenX,lenY))
spline = biCubicSpline(xdata, ydata, zdata, 20, 0.0001)

x = range(xdata[1], xdata[end], N)
y = range(ydata[1], ydata[end], M)
z = evaluate(spline, N, M)

# Create the plot
#gr()
xdata = xPlotData_1
ydata = yPlotData_1
zdata = zPlotData_1

plot(x, y, z, seriestype=:surface, label="Lavery spline")
plot!(xdata, ydata, zdata, seriestype=:scatter, label="data points", camera=(-90,90))
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

plot(x, y, z, seriestype=:surface, label="Lavery spline")
plot!(xdata, ydata, zdata, seriestype=:scatter, label="data points", camera=(90,0))
xaxis!("x")  # Set x-axis label
yaxis!("y")  # Set y-axis label
zaxis!("z")  # Set z-axis label
title!("Points interpolated with Bi-variate Lavery Splines")
#display(p)
savefig("bi_lavery_side.png")