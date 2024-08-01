using Plots
using DelimitedFiles

include("laverySplinesBi.jl")
include("testData.jl")
include("plotData.jl")

xdata = xdata_5
ydata = ydata_5
zdata = zdata_5
lenX = length(xdata)
lenY = length(ydata)

# spline = Data(xdata, ydata, zdata, zeros(lenX,lenY), zeros(lenX,lenY))
spline = biCubicSpline(xdata, ydata, zdata, 20, 0.0001)

file = open("bx_data.csv", "w")
writedlm(file, spline.bx, ',')
close(file)
file = open("by_data.csv", "w")
writedlm(file, spline.by,',')
close(file)

N = 30
M = N
x = range(xdata[1], xdata[end], N)
y = range(ydata[1], ydata[end], M)
z = evaluate(spline, N, M)

p = plot3D(xdata, ydata, zdata,z,save_gif=true)

savefig(p,"3Dplot.png")

