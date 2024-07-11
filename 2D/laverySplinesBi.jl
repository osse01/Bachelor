
using JuMP, HiGHS
model = Model(HiGHS.Optimizer)

function biCubicSpline(xData, yData, zData)
    len = length(xData)                             # length of xData
    deltaX = xData[i+1]-xData[i]                    # x-step length
    deltaY = yData[i+1]-yData[i]                    # y-step length

    @variable(model, dzdx[1:len,1:len])             # dzdx(x_i,y_i)
    @variable(model, dzdy[1:len,1:len])             # dzdy(x_i,y_i)

    d2zdx2(i,j,k,l,xTilde, yTilde) = 1/(deltaX^2)*((-6+12*xTilde)*zData[i,j]
                                    +deltaX*(-4+6*xTilde)*dzdx[i,j]
                                    +(6-12*xTilde)*zData[i+1,j]
                                    +deltaX*(-1+6*deltaX)*dzdx[i+1,j])
    d2zdxdy(i,j,k,l,xTilde, yTilde) = 1/(deltaX*deltaY)*(6*yTilde*(zData[i,j]-zData[i+1,j])
                                    +deltaX*yTilde*dzdx[i,j]
                                    +deltaY*(-1+2*yTilde)*dzdy[i,j]
                                    
                                        ) 

    @objective(model, Min, sum(sum(1/(N^2)(
        sum(sum(gamma(arg1)+gamma(arg2)+gamma(arg3)+gamma(arg4)))
        +sum(sum(gamma(arg1)+gamma(arg2)+gamma(arg3)+gamma(arg4)))))))
    return
end