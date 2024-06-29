using Plots
#using CubicSplines
include("laverySpline.jl")
#using Interpolations

len=50 # number of interpolations

# Generate sinusoidal data with randomness in the x coordinates
xdata = range(0, stop=15, length=len) .+ 10/len*rand(len)
ydata = sin.(xdata)
#ydata = vec( [zeros( 1, Int(len/2) ) ones( 1, Int(len/2) )] )
#ydata = vec( ones(1, Int(len)) )

# Evaluate the spline at more points for smoother plotting
xs = range(xdata[1], stop=xdata[end], length=1000)
intervall = vec(xs)
ys = laverySpline(xdata, ydata, intervall)
print(ys)
#derivative = [gradient(spline, x) for x in xs] # First derivative

# Create the plot
plot(xdata, ydata, seriestype=:scatter, label="Data points")
plot!(xs, ys, label="Lavery spline")
#plot!(xs, derivative, label="Derivative")
xaxis!("x")  # Set x-axis label
yaxis!("y")  # Set y-axis label
title!("Function with Lavery Splines")
#legend!()
ylims = (0,2)

savefig("data.png")
