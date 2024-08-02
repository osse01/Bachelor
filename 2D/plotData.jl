

# Create the plot
function plot3D(xdata,ydata,zdata,z; show_data=true, camera_angle=(30,30), title::AbstractString="Beautiful Plots!")
    if show_data
        x_scatter_data = repeat(xdata,length(ydata))
        y_scatter_data = repeat(ydata,inner=length(xdata))
        z_scatter_data = vec(zdata')
    end
    p = plot(x, y, z, seriestype=:surface, label="Lavery spline",camera=camera_angle)
    if show_data
        p = plot!(x_scatter_data, y_scatter_data, z_scatter_data, seriestype=:scatter, label="data points")
    end
    xaxis!("x")  # Set x-axis label
    yaxis!("y")  # Set y-axis label
    zaxis!("z")  # Set z-axis label
    title!(title)
    return p
end
function saveGIF(xdata,ydata,zdata,z; show_data=true, title::AbstractString="GIFs, more like GILFS")

    if show_data
        x_scatter_data = repeat(xdata,length(ydata))
        y_scatter_data = repeat(ydata,inner=length(xdata))
        z_scatter_data = vec(zdata')
    end
    gr()
    anim = Animation()
    for i in 0:359
        g = plot(x, y, z, seriestype=:surface, label="Lavery spline",camera=(i,30))
        if show_data
            g = plot!(g, x_scatter_data, y_scatter_data, z_scatter_data, seriestype=:scatter, label="data points",camera=(i,30))
        end
        xaxis!("x")  # Set x-axis label
        yaxis!("y")  # Set y-axis label
        zaxis!("z")  # Set z-axis label
        title!(title)
        frame(anim,g) 
    end
    
    return anim
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