

# Create the plot
function plot3D(xdata,ydata,zdata,z, show_data=true)
    N = 100
    M = N
    x = range(xdata[1], xdata[end], N)
    y = range(ydata[1], ydata[end], M)
    x_scatter_data = repeat(xdata,length(ydata))
    y_scatter_data = repeat(ydata,inner=length(xdata))
    z_scatter_data = vec(zdata')
    p = plot(x, y, z, seriestype=:surface, label="Lavery spline")
    if show_data
        p = plot!(x_scatter_data, y_scatter_data, z_scatter_data, seriestype=:scatter, label="data points")
    end
    xaxis!("x")  # Set x-axis label
    yaxis!("y")  # Set y-axis label
    zaxis!("z")  # Set z-axis label
    title!("Points interpolated with Bi-variate Lavery Splines")
    return p
end

# function plot2()
#     return
# end

# function plot3()
#     plot(x, y, z, seriestype=:surface, label="Lavery spline")
#     plot!(xdata, ydata, zdata, seriestype=:scatter, label="data points", camera=(-90,90))
#     xaxis!("x")  # Set x-axis label
#     yaxis!("y")  # Set y-axis label
#     zaxis!("z")  # Set z-axis label
#     title!("Points interpolated with Bi-variate Lavery Splines")

#     p  = plot(x, y, z, seriestype=:surface, label="Lavery spline")
#     plot!(xdata, ydata, zdata, seriestype=:scatter, label="data points", camera=(90,0))
#     xaxis!("x")  # Set x-axis label
#     yaxis!("y")  # Set y-axis label
#     zaxis!("z")  # Set z-axis label
#     title!("Points interpolated with Bi-variate Lavery Splines")

#     return p
# end