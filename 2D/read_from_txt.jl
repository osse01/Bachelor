
struct txtData
    x::Vector
    y::Vector
    z::Array
end


function read_values_from_file(file_path::AbstractString)
        # Open the file
        file = open(file_path, "r")

        # Initialize an empty array to store the values
        values = Vector{Float64}()

        # Read the remaining lines
        for line in eachline(file)
            # Ignore lines starting with #
            if !startswith(line, "#")
                # Convert the line to a float and append it to the array
                push!(values, parse(Float64, line))
            end
        end

        # Close the file
        close(file)

        # Extract number of x-values and y-values
        n = Integer(popfirst!(values))
        m = Integer(popfirst!(values))
        
        # Extract the x-values and y-values
        x_vals = values[1:n]
        y_vals = values[n+1:n+m]
        z_vals = values[n+1+m:end]
        print(length(x_vals))
        print(length(y_vals))
        print(length(z_vals))

        z_array = reshape(z_vals,n,m)

        data = txtData(x_vals,y_vals,z_array)

        # Return the fdata
        return data
end