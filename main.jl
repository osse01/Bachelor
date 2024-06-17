using Plots
using CubicSplines
#using Interpolations

len=50 # number of interpolations

# Generate sinusoidal data with randomness in the x coordinates
xdata = range(0, stop=15, length=len) .+ 10/len*rand(len)
#ydata = sin.(xdata) .+ 0.1rand(20)
ydata = [zeros(Int(len/2),1);
         ones(Int(len/2),1)]

# Create a CubicSpline object
spline = CubicSpline(xdata, ydata)


# Evaluate the spline at more points for smoother plotting
xs = range(xdata[1], stop=xdata[end], length=1000)
ys = spline(xs)
derivative = [gradient(spline, x) for x in xs] # First derivative

# Create the plot
plot(xdata, ydata, seriestype=:scatter, label="Data points")
plot!(xs, ys, label="Cubic spline")
plot!(xs, derivative, label="Derivative")
xaxis!("x")  # Set x-axis label
yaxis!("y")  # Set y-axis label
title!("Sinusoidal Data with Cubic Spline")
#legend!()

savefig("data.png")
