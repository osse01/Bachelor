
using JuMP, HiGHS
model = Model(HiGHS.Optimizer)

# Calculates the Lavery Splines of xdata and zdata
function laverySpline(xData,zData, intervall)
    len = length(xData)
    h(i) = xData[i+1]-xData[i]
    deltaZ(i) = (zData[i+1]-zData[i]) / h(i)
    weight = 1/len
    
    # HiGHS solver using JuMP
    @variable(model, b[1:len])
    bDomain = range(1,len)
    integralSteps = Int(300)
    integralDomain = range(-0.5, 0.5, integralSteps)
    sumDomain = range(1,len - 1)
    @variable(model, abs_E[sumDomain, integralDomain]>=0) # Absolute value
    @variable(model, abs_b[1:len])
    integralArgument(i,t) = (-1+6*t)*b[i] + (1+6*t)*b[i+1] - 12*deltaZ(i)*t

    @objective(model, Min, sum( sum( abs_E[i,t] / integralSteps[i]
                                for t in integralDomain)
                                    for i in  sumDomain) +
                             sum( weight * abs_b[i] for i in bDomain))
    @constraint(model, argPos[t in integralDomain, i in sumDomain],
                                        abs_E[i, t] >= integralArgument(i, t) )
    @constraint(model, argNeg[t in integralDomain, i in sumDomain],
                                        abs_E[i, t] >= -integralArgument(i, t) )
    @constraint(model, bPos[i in bDomain], abs_b[i] >= b[i])
    @constraint(model, bNeg[i in bDomain], abs_b[i] >= -b[i])
    
    optimize!(model)
    b = vec(value.(b))

    x(j) = intervall[j]
    zTerm(i,j) = zData[i] + b[i]*(x(j) - xData[i]) + 
                1/h(i) * (-(2*b[i]+b[i+1]) + 3*deltaZ(i))*(x(j)-xData[i])^2 + 
                1/h(i)^2 * (b[i]+b[i+1]-2*deltaZ(i))*(x(j)-xData[i])^3


    z = zeros(length(intervall),1)
    i = 1
    for j in range(1,length(intervall))
        if x(j) > xData[i+1]
            i = i+1
        end
    z[j] = zTerm(i,j)
    end
    return z
end