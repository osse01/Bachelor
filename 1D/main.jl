using Plots
#using Interpolations
using BasicInterpolators

include("laverySpline.jl")
include("testData.jl")

#len=60 # number of interpolations

# Generate data with randomness in the x coordinates:
#xdata = range(0, stop=15, length=len) .+ 10/len*rand(len)
# Set ydata:
#ydata = sin.(xdata)
#ydata = vec( [zeros( 1, Int(len/2) ) 100*ones( 1, Int(len/2) )] )
#ydata = vec( ones(1, Int(len)) )
#ydata = vec( [zeros( 1, Int(len/3) ) ones( 1, Int(len/3) ) zeros( 1, Int(len/3) )] )


xdata = xDataRect
ydata = yDataRect

# Evaluate the spline at more points for smoother plotting
xs = range(xdata[1], stop=xdata[end], length=1000)
ys = laverySpline(xdata, ydata, xs)

# Create the plot
plot(xdata, ydata, seriestype=:scatter, label="Data points")
plot!(xs, ys, label="Lavery spline")
xaxis!("x")  # Set x-axis label
yaxis!("y")  # Set y-axis label
title!("Function with Lavery Splines")
savefig("Lavery_rect.png")

#spline = CubicSpline(xDataMedium, yDataMedium)
itp  = CubicSplineInterpolator(xdata, ydata, WeakBoundaries());
#interpolate((xdata,), ydata, Gridded(Cubic(Line(OnGrid()))))
xs = range(xdata[1], stop=xdata[end], length=1000)
ys = [itp(x) for x in xs]
plot(xdata, ydata, seriestype=:scatter, label="Data points")
plot!(xs, ys, label="Cubic spline")
xaxis!("x")  # Set x-axis label
yaxis!("y")  # Set y-axis label
title!("Function with Cubic Splines")
savefig("Cubic_rect.png")
