using Plots
include("laverySplinesBi.jl")
include("testData.jl")
include("plotData.jl")

xdata = xdata_4
ydata = ydata_4
zdata = zdata_4
lenX = length(xdata)
lenY = length(ydata)

N = 100
M = N
#spline = Data(xdata, ydata, zdata, zeros(lenX,lenY), zeros(lenX,lenY))
spline = biCubicSpline(xdata, ydata, zdata, 20, 0.0001)

file = open("bx_data.csv", "w")
write(file, splines.bx)
close(file)
file = open("by_data.csv", "w")
write(file, splines.by)
close(file)

z = evaluate(spline, N, M)

p = plot3D(xdata, ydata, zdata,z,true)

savefig(p,"3Dplot.png")

