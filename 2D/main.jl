using Plots
using DelimitedFiles
using BasicInterpolators

include("laverySplinesBi.jl")
include("testData.jl")
include("plotData.jl")


xdata = xdata_7
ydata = ydata_7
zdata = zdata_7
lenX = length(xdata)
lenY = length(ydata)


#bx = Array(readdlm("bx_data_6.csv", ',', Float64))
#by = Array(readdlm("by_data_6.csv", ',', Float64))

#bx = zeros(lenX,lenY)
#by = zeros(lenX,lenY)
#spline = Data(xdata, ydata, zdata, bx, by)
print("\nStart Interpolating...\n")
spline, model = biCubicSpline(xdata, ydata, zdata, 1, 0)

#print(model)
#file = open("modelOutput", "w")
#write(file, model)
#close(file)

# Saving derivatives to file
print("Print bx. by to file...\n")
file = open("bx_data_7.csv", "w")
writedlm(file, spline.bx, ',')
close(file)
file = open("by_data_7.csv", "w")
writedlm(file, spline.by,',')
close(file)

N = 500
M = N
x = range(xdata[1], xdata[end], N)
y = range(ydata[1], ydata[end], M)
z = evaluate(spline, N, M)
print("Plotting...\n")
p = plot3D(xdata, ydata, zdata,x,y,z,show_data=false,title="Lavery Splines")
savefig(p,"geo_Lavery.png")

anim = saveGIF(xdata, ydata, zdata,z,show_data=false,title="Lavery Splines")
gif(anim,"geo_Lavery.gif",fps=30)

#BiCubic
spline = BicubicInterpolator(xdata, ydata, zdata)
xs = range(xdata[1], xdata[end], N)
ys = range(ydata[1], ydata[end], M)
z_cubic = [spline(x, y) for x in xs, y in ys]
print(size(z))
p = plot3D(xdata, ydata, zdata,x,y,z_cubic,show_data=false, title="Bi-Cubic Splines")
savefig(p,"geo_Bicubic.png")

p = plot3D(xdata, ydata, zdata,x,y,z_cubic-z,show_data=false, title="Diffrence in Splines")
savefig(p,"geo_diff.png")

anim = saveGIF(xdata, ydata, zdata,z,show_data=false,title="Cubic Splines")
gif(anim,"geo_Bicubic.gif",fps=30)

