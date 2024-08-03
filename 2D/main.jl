using Plots
using DelimitedFiles
using BasicInterpolators

include("laverySplinesBi.jl")
include("testData.jl")
include("plotData.jl")


xdata = xdata_6
ydata = ydata_6
zdata = zdata_6
lenX = length(xdata)
lenY = length(ydata)


bx = Array(readdlm("bx_data_6.csv", ',', Float64))
by = Array(readdlm("by_data_6.csv", ',', Float64))

#bx = zeros(lenX,lenY)
#by = zeros(lenX,lenY)
spline = Data(xdata, ydata, zdata, bx, by)
print("\nStart Interpolating...\n")
#spline = biCubicSpline(xdata, ydata, zdata, 5, 0.0001)

# file = open("bx_data.csv", "w")
# writedlm(file, spline.bx, ',')
# close(file)
# file = open("by_data.csv", "w")
# writedlm(file, spline.by,',')
# close(file)

N = 1000
M = N
x = range(xdata[1], xdata[end], N)
y = range(ydata[1], ydata[end], M)
z = evaluate(spline, N, M)

p = plot3D(xdata, ydata, zdata,z,show_data=false,title="LaverySplines are nice")
savefig(p,"3Dplot.png")

anim = saveGIF(xdata, ydata, zdata,z,show_data=false,title="GIFs are cool")
gif(anim,"spline.gif",fps=30)

#BiCubic
#spline = BicubicInterpolator(xdata, ydata, zdata)
#xs = range(xdata[1], xdata[end], N)
#ys = range(ydata[1], ydata[end], M)
#z = [spline(x, y) for x in xs, y in ys]
#print(size(z))
#p = plot3D(xdata, ydata, zdata,z,show_data=false, save_gif=true)

#savefig(p,"3Dplot_Bicubic.png")

