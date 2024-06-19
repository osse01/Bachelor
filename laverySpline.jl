
using JuMP
using HiGHS
using Integrals

# Evaluate E given parameters recursively
function E(x,z,b)
    if (length(x) == 2) 
        return 0
    else
    deltaZ = (z[2] - z[1]) / (x[2] - x[1])
    domain = (-0.5, 0.5)
    func(t, p) = abs( (-1+6*t)*b[1] + (1+6*t)*b[2] - 12*deltaZ*t )
    p=1
    popfirst!(x)
    popfirst!(z)
    popfirst!(b)
    prob = IntegralProblem(func, domain)
    sol = solve(prob, QuadGKJL(); reltol = 1e-3)
    return E(x,z,b) + sol.u
    end
end

# Calculates the Lavery Splines of xdata and zdata
function laverySpline(xEval,xData,zData)
    len = length(xData)
    h(i) = xData[i+1]-xData[i]
    deltaZ(i) = (zData[i+1]-zData[i]) / h(i)
    b = vec([zeros(Int(len/2),1); ones(Int(len/2),1)])
    
    #E(xData,zData,b) # Use to minimize the derivative oscillations
    
    z(i) = zData[i] + b[i]*(xEval - xData[i]) + 1/h(i) * (-(2*b[i]+b[i+1]) + 3*deltaZ(i))*(xEval-xData[i])^2 + 1/h(i)^2 * (b[i]+b[i+1]-2*deltaZ(i))*(xEval-xData[i])^3

    # Which index?
    for i in 1:(length(xData)-1)
        if xData[i] <= xEval && xEval <= xData[i+1]
            return z(i)
        end
    end
end