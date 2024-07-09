using Plots
include("laverySpline.jl")

len=9 # number of interpolations

# Generate data with randomness in the x coordinates:
xdata = range(0, stop=15, length=len) .+ 10/len*rand(len)
# Set ydata:
ydata = sin.(xdata)
#ydata = vec( [zeros( 1, Int(len/2) ) 100*ones( 1, Int(len/2) )] )
#ydata = vec( ones(1, Int(len)) )
#ydata = vec( [zeros( 1, Int(len/3) ) ones( 1, Int(len/3) ) zeros( 1, Int(len/3) )] )


# Evaluate the spline at more points for smoother plotting
xs = range(xdata[1], stop=xdata[end], length=1000)
intervall = vec(xs)
ys = laverySpline(xdata, ydata, intervall)
print(ys)

# Create the plot
plot(xdata, ydata, seriestype=:scatter, label="Data points")
plot!(xs, ys, label="Lavery spline")
xaxis!("x")  # Set x-axis label
yaxis!("y")  # Set y-axis label
title!("Function with Lavery Splines")

savefig("data.png")
