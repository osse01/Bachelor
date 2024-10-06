
using Plots
include("testData.jl")

xData = xDataStep
zData = yDataStep

# Evaluate the spline at more points for smoother plotting
xs = range(xData[1], stop=xData[end], length=1000)
intervall = xs

    len = length(xData)
    h(i) = xData[i+1]-xData[i]
    deltaZ(i) = (zData[i+1]-zData[i]) / h(i)

    b = vec([zeros(10,1)])

    x(j) = intervall[j]
    zTerm(i,j) = zData[i] + b[i]*(x(j) - xData[i]) + 
                1/h(i) * (-(2*b[i]+b[i+1]) + 3*deltaZ(i))*(x(j)-xData[i])^2 + 
                1/h(i)^2 * (b[i]+b[i+1]-2*deltaZ(i))*(x(j)-xData[i])^3


    z = zeros(length(intervall),1)
    global k = 1
    for j in range(1,length(intervall))
        if x(j) > xData[k+1]
            global k = k+1
        end
        z[j] = zTerm(k,j)
    end




# Create the plot
plot(xdata, ydata, seriestype=:scatter, label="Data points")
plot!(xs, ys, label="Lavery spline")
xaxis!("x")  # Set x-axis label
yaxis!("y")  # Set y-axis label
title!("Function with Lavery Splines")
savefig("Small_step.png")