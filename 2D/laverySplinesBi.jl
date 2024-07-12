
using JuMP, HiGHS
model = Model(HiGHS.Optimizer)

# set_attribute(model, "presolve", "on")
# set_attribute(model, "solver", "ipm") # Interior Point Method

function biCubicSpline(xData, yData, zData, N)
    I = length(xData)                               # length of xData
    J = length(yData)                               # length of yData
    deltaX = [xData[i+1]-xData[i] for i in 1:I-1]   # x-step length
    deltaY = [yData[j+1]-yData[j] for j in 1:J-1]   # y-step length
    xTilde(x,i) = (x - xData[i]) / deltaX[i]        # Transforms to unit cube
    yTilde(y,i) = (y - yData[i]) / deltaY(i)        # Transforms to unit cube

    @variable(model, bx[1:len,1:len])               # dzdx(x_i,y_i)
    @variable(model, by[1:len,1:len])               # dzdy(x_i,y_i)

    d2zdx2_1(i,j,k,l,xTilde, yTilde) = 1/(deltaX[i]^2)*((-6+12*xTilde)*zData[i,j]
                                    +deltaX[i]*(-4+6*xTilde)*bx[i,j]
                                    +(6-12*xTilde)*zData[i+1,j]
                                    +deltaX[i]*(-1+6*deltaX[i])*bx[i+1,j])
    d2z1dxdy_1(i,j,k,l,xTilde, yTilde) = 1/(deltaX[i]*deltaY[j])*(
                                    6*yTilde*( (zData[i,j]+zData[i+1,j+1]) - (zData[i+1,j]+zData[i,j+1]) )
                                    +deltaX[i]*yTilde*( (bx[i,j]+bx[i+1,j]) - (bx[i,j+1]+bx[i+1,j+1]) )
                                    +deltaY[j]*( (by[i1,j]-by[i,j]) + 
                                    2*yTilde*( (by[i,j]+by[i,j+1]) - (by[i+1,j]+by[i+1,j+1]) ) ) )
    d2z1dy2_1(i,j,k,l,xTilde,yTilde) = 1/(deltaX[i]^2)*((-6+6*xTilde+6*yTilde)*zData[i,j]
                                    +deltaX[i]*(-1+xTilde)*bx[i,j]
                                    +deltaY[j]*(-3+2*xTilde+3*yTilde)*by[i,j]
                                    +(-6*xTilde+6*yTilde)*zData[i+1,j]
                                    +deltaX[i]*(xTilde)*bx[i+1,j]
                                    +deltaY[j]*(-1-2*xTilde+3*yTilde)*by[i+1,j]
                                    +(6-6*xTilde-6*yTilde)*zData[i,j+1]
                                    +deltaX[i]*(1-xTilde)*bx[i,j+1]
                                    +deltaY[j]*(-2+2*xTilde+3*yTilde)*by[i,j+1]
                                    +(6*xTilde-6*yTilde)*zData[i+1,j+1]
                                    +deltaX[i]*(-xTilde)*bx[i+1,j+1]
                                    +deltaY[j]*(-2*xTilde+3*yTilde)by[i+1,j+1]
                                    )

    d2zdx2_2(i,j,k,l,xTilde, yTilde) = 1/(deltaX[i]^2)*((-6+12*xTilde)*zData[i,j]
                                    +deltaX[i]*(-4+6*xTilde)*bx[i,j]
                                    +(6-12*xTilde)*zData[i+1,j]
                                    +deltaX[i]*(-1+6*deltaX[i])*bx[i+1,j])
    d2z1dxdy_2(i,j,k,l,xTilde, yTilde) = 1/(deltaX[i]*deltaY[j])*(
                                    6*yTilde*( (zData[i,j]+zData[i+1,j+1]) - (zData[i+1,j]+zData[i,j+1]) )
                                    +deltaX[i]*yTilde*( (bx[i,j]+bx[i+1,j]) - (bx[i,j+1]+bx[i+1,j+1]) )
                                    +deltaY[j]*( (by[i1,j]-by[i,j]) + 
                                    2*yTilde*( (by[i,j]+by[i,j+1]) - (by[i+1,j]+by[i+1,j+1]) ) ) )
    d2z1dy2_2(i,j,k,l,xTilde,yTilde) = 1/(deltaX[i]^2)*((-6+6*xTilde+6*yTilde)*zData[i,j]
                                    +deltaX[i]*(-1+xTilde)*bx[i,j]
                                    +deltaY[j]*(-3+2*xTilde+3*yTilde)*by[i,j]
                                    +(-6*xTilde+6*yTilde)*zData[i+1,j]
                                    +deltaX[i]*(xTilde)*bx[i+1,j]
                                    +deltaY[j]*(-1-2*xTilde+3*yTilde)*by[i+1,j]
                                    +(6-6*xTilde-6*yTilde)*zData[i,j+1]
                                    +deltaX[i]*(1-xTilde)*bx[i,j+1]
                                    +deltaY[j]*(-2+2*xTilde+3*yTilde)*by[i,j+1]
                                    +(6*xTilde-6*yTilde)*zData[i+1,j+1]
                                    +deltaX[i]*(-xTilde)*bx[i+1,j+1]
                                    +deltaY[j]*(-2*xTilde+3*yTilde)by[i+1,j+1]
                                    )
    d2zdx2_3(i,j,k,l,xTilde, yTilde) = 1/(deltaX[i]^2)*((-6+12*xTilde)*zData[i,j]
                                    +deltaX[i]*(-4+6*xTilde)*bx[i,j]
                                    +(6-12*xTilde)*zData[i+1,j]
                                    +deltaX[i]*(-1+6*deltaX[i])*bx[i+1,j])
    d2z1dxdy_3(i,j,k,l,xTilde, yTilde) = 1/(deltaX[i]*deltaY[j])*(
                                    6*yTilde*( (zData[i,j]+zData[i+1,j+1]) - (zData[i+1,j]+zData[i,j+1]) )
                                    +deltaX[i]*yTilde*( (bx[i,j]+bx[i+1,j]) - (bx[i,j+1]+bx[i+1,j+1]) )
                                    +deltaY[j]*( (by[i1,j]-by[i,j]) + 
                                    2*yTilde*( (by[i,j]+by[i,j+1]) - (by[i+1,j]+by[i+1,j+1]) ) ) )
    d2z1dy2_3(i,j,k,l,xTilde,yTilde) = 1/(deltaX[i]^2)*((-6+6*xTilde+6*yTilde)*zData[i,j]
                                    +deltaX[i]*(-1+xTilde)*bx[i,j]
                                    +deltaY[j]*(-3+2*xTilde+3*yTilde)*by[i,j]
                                    +(-6*xTilde+6*yTilde)*zData[i+1,j]
                                    +deltaX[i]*(xTilde)*bx[i+1,j]
                                    +deltaY[j]*(-1-2*xTilde+3*yTilde)*by[i+1,j]
                                    +(6-6*xTilde-6*yTilde)*zData[i,j+1]
                                    +deltaX[i]*(1-xTilde)*bx[i,j+1]
                                    +deltaY[j]*(-2+2*xTilde+3*yTilde)*by[i,j+1]
                                    +(6*xTilde-6*yTilde)*zData[i+1,j+1]
                                    +deltaX[i]*(-xTilde)*bx[i+1,j+1]
                                    +deltaY[j]*(-2*xTilde+3*yTilde)by[i+1,j+1]
                                    )
    d2zdx2_4(i,j,k,l,xTilde, yTilde) = 1/(deltaX[i]^2)*((-6+12*xTilde)*zData[i,j]
                                    +deltaX[i]*(-4+6*xTilde)*bx[i,j]
                                    +(6-12*xTilde)*zData[i+1,j]
                                    +deltaX[i]*(-1+6*deltaX[i])*bx[i+1,j])
    d2z1dxdy_4(i,j,k,l,xTilde, yTilde) = 1/(deltaX[i]*deltaY[j])*(
                                    6*yTilde*( (zData[i,j]+zData[i+1,j+1]) - (zData[i+1,j]+zData[i,j+1]) )
                                    +deltaX[i]*yTilde*( (bx[i,j]+bx[i+1,j]) - (bx[i,j+1]+bx[i+1,j+1]) )
                                    +deltaY[j]*( (by[i1,j]-by[i,j]) + 
                                    2*yTilde*( (by[i,j]+by[i,j+1]) - (by[i+1,j]+by[i+1,j+1]) ) ) )
    d2z1dy2_4(i,j,k,l,xTilde,yTilde) = 1/(deltaX[i]^2)*((-6+6*xTilde+6*yTilde)*zData[i,j]
                                    +deltaX[i]*(-1+xTilde)*bx[i,j]
                                    +deltaY[j]*(-3+2*xTilde+3*yTilde)*by[i,j]
                                    +(-6*xTilde+6*yTilde)*zData[i+1,j]
                                    +deltaX[i]*(xTilde)*bx[i+1,j]
                                    +deltaY[j]*(-1-2*xTilde+3*yTilde)*by[i+1,j]
                                    +(6-6*xTilde-6*yTilde)*zData[i,j+1]
                                    +deltaX[i]*(1-xTilde)*bx[i,j+1]
                                    +deltaY[j]*(-2+2*xTilde+3*yTilde)*by[i,j+1]
                                    +(6*xTilde-6*yTilde)*zData[i+1,j+1]
                                    +deltaX[i]*(-xTilde)*bx[i+1,j+1]
                                    +deltaY[j]*(-2*xTilde+3*yTilde)by[i+1,j+1]
                                    )

    z1(x,y,bx,by,i,j) =
        (1 - 3*xTilde(x,i)^2 + 2*xTilde(x,i) - 3*yTilde(y,j) + 3*xTilde(x,i)*yTilde(y,j)^2 + yTilde(y,j)^3)*zData[i,j] +
        deltaX[i] * (xTilde(x,i) - 2*xTilde(x,i)^2 + xTilde(x,i)^3 - 0.5 * xTilde(x,i)*yTilde(y,j)^2) * bx[i,j] + 
        deltaY[j] * (yTilde(y,j) - xTilde(x,i)*yTilde(y,j) - 1.5 * yTilde(y,j)^2 + xTilde(x,i)*yTilde(y,j)^2 + 0.5*yTilde(y,j)^3) * by[i,j] +
        (3*xTilde(x,i)^2 - 2xTilde(x,i)^3 - 3*xTilde(x,i)*yTilde(y,j)^2 + yTilde(y,j)^3) * zData[i+1,j] +
        deltaX[i] * (-xTilde(x,i)^2 + xTilde(x,i)^3 + 0.5*xTilde(x,i)* yTilde(y,j)^2) * bx[i+1,j] +
        deltaY[j] * (xTilde(x,i) * yTilde(y,j) - 0.5*yTilde(y,j)^2 - xTilde(x,i)*yTilde(y,j)^2 + 0.5*yTilde(y,j)^3) * by[i+1,j] +
        (3*yTilde(y,j)^2 - 3*xTilde(x,i) * yTilde(y,j)^2 - yTilde(y,j)^3)*zData[i,j+1] + deltaX[i]*(0.5*yTilde(y,j)^2-0.5*xTilde(x,i)*yTilde(y,j))*bx[i,j+1] +
        deltaY[j] * (-yTilde(y,j)^2 + xTilde(x,i)*yTilde(y,j)^2 + 0.5*yTilde(y,j)^3)*by[i,j+1] + (3*xTilde(x,i)*yTilde(y,j)^2 - yTilde(y,j)^3)*zData[i+1,j+1] +
        deltaX[i] * (-0.5xTilde(x,i)*yTilde(y,j)^2)* bx[i+1,j+1] + deltaY[j]*(-xTilde(x,i)*yTilde(y,j)^2 + 0.5*yTilde(y,j)^3)*by[i+1,j+1]

    z2(x,y,bx,by,i,j) =
        (1 - 3*yTilde(y,j)^2 + 2*yTilde(y,j) - 3*( 1-xTilde(x,i) ) + 3*yTilde(y,j)*( 1-xTilde(x,i) )^2 + ( 1-xTilde(x,i) )^3)*zData[i+1,j] +
        deltaY[j] * (yTilde(y,j) - 2*yTilde(y,j)^2 + yTilde(y,j)^3 - 0.5 * yTilde(y,j)*( 1-xTilde(x,i) )^2) * by[i+1,j] + 
        -deltaX[i] * (( 1-xTilde(x,i) ) - yTilde(y,j)*( 1-xTilde(x,i) ) - 1.5 * ( 1-xTilde(x,i) )^2 + yTilde(y,j)*( 1-xTilde(x,i) )^2 + 0.5*( 1-xTilde(x,i) )^3) * -bx[i+1,j] +
        (3*yTilde(y,j)^2 - 2yTilde(y,j)^3 - 3*yTilde(y,j)*( 1-xTilde(x,i) )^2 + ( 1-xTilde(x,i) )^3) * zData[i+1,j+1] +
        deltaY[j] * (-yTilde(y,j)^2 + yTilde(y,j)^3 + 0.5*yTilde(y,j)* ( 1-xTilde(x,i) )^2) * by[i+1,j+1] +
        -deltaX[i] * (yTilde(y,j) * ( 1-xTilde(x,i) ) - 0.5*( 1-xTilde(x,i) )^2 - yTilde(y,j)*( 1-xTilde(x,i) )^2 + 0.5*( 1-xTilde(x,i) )^3) * -bx[i+1,j+1] +
        (3*( 1-xTilde(x,i) )^2 - 3*yTilde(y,j) * ( 1-xTilde(x,i) )^2 - ( 1-xTilde(x,i) )^3)*zData[i,j] + deltaY[j]*(0.5*( 1-xTilde(x,i) )^2-0.5*yTilde(y,j)*( 1-xTilde(x,i) ))*by[i,j] +
        -deltaX[i] * (-( 1-xTilde(x,i) )^2 + yTilde(y,j)*( 1-xTilde(x,i) )^2 + 0.5*( 1-xTilde(x,i) )^3)*-bx[i,j] + (3*yTilde(y,j)*( 1-xTilde(x,i) )^2 - ( 1-xTilde(x,i) )^3)*zData[i,j+1] +
        deltaY[j] * (-0.5*yTilde(y,j)*( 1-xTilde(x,i) )^2)* by[i,j+1] - deltaX[i]*(-yTilde(y,j)*( 1-xTilde(x,i) )^2 + 0.5*( 1-xTilde(x,i) )^3)*-bx[i,j+1]
    
    z3(x,y,bx,by,i,j) =
        (1 - 3*( 1-xTilde(x,i) )^2 + 2*( 1-xTilde(x,i) ) - 3*( 1-yTilde(y,j) ) + 3*( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2 + ( 1-yTilde(y,j) )^3)*zData[i+1,j+1] +
        -deltaX[i] * (( 1-xTilde(x,i) ) - 2*( 1-xTilde(x,i) )^2 + ( 1-xTilde(x,i) )^3 - 0.5 * ( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2) * -bx[i+1,j+1] + 
        -deltaY[j] * (( 1-yTilde(y,j) ) - ( 1-xTilde(x,i) )*( 1-yTilde(y,j) ) - 1.5 * ( 1-yTilde(y,j) )^2 + ( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2 + 0.5*( 1-yTilde(y,j) )^3) * -by[i+1,j+1] +
        (3*( 1-xTilde(x,i) )^2 - 2( 1-xTilde(x,i) )^3 - 3*( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2 + ( 1-yTilde(y,j) )^3) * zData[i,j+1] +
        -deltaX[i] * (-( 1-xTilde(x,i) )^2 + ( 1-xTilde(x,i) )^3 + 0.5*( 1-xTilde(x,i) )* ( 1-yTilde(y,j) )^2) * -bx[i,j+1] +
        -deltaY[j] * (( 1-xTilde(x,i) ) * ( 1-yTilde(y,j) ) - 0.5*( 1-yTilde(y,j) )^2 - ( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2 + 0.5*( 1-yTilde(y,j) )^3) * -by[i,j+1] +
        (3*( 1-yTilde(y,j) )^2 - 3*( 1-xTilde(x,i) ) * ( 1-yTilde(y,j) )^2 - ( 1-yTilde(y,j) )^3)*zData[i+1,j] + -deltaX[i]*(0.5*( 1-yTilde(y,j) )^2-0.5*( 1-xTilde(x,i) )*( 1-yTilde(y,j) ))*-bx[i+1,j] +
        -deltaY[j] * (-( 1-yTilde(y,j) )^2 + ( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2 + 0.5*( 1-yTilde(y,j) )^3)*-by[i+1,j] + (3*( 1-xTilde(x,i) )*( -yTilde(y,j) )^2 - ( 1-yTilde(y,j) )^3)*zData[i,j] +
        -deltaX[i] * (-0.5( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2)* -bx[i,j] + -deltaY[j]*(-( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2 + 0.5*( 1-yTilde(y,j) )^3)*-by[i,j]

    z4(x,y,bx,by,i,j) =
        (1 - 3*( 1-yTilde(y,j) )^2 + 2*( 1-yTilde(y,j) ) - 3*xTilde(x,i) + 3*( 1-yTilde(y,j) )*xTilde(x,i)^2 + xTilde(x,i)^3)*zData[i,j+1] +
        -deltaY[j] * (( 1-yTilde(y,j) ) - 2*( 1-yTilde(y,j) )^2 + ( 1-yTilde(y,j) )^3 - 0.5 * ( 1-yTilde(y,j) )*xTilde(x,i)^2) * -by[i,j+1] + 
        deltaX[i] * (xTilde(x,i) - ( 1-yTilde(y,j) )*xTilde(x,i) - 1.5 * xTilde(x,i)^2 + ( 1-yTilde(y,j) )*xTilde(x,i)^2 + 0.5*xTilde(x,i)^3) * bx[i,j+1] +
        (3*( 1-yTilde(y,j) )^2 - 2( 1-yTilde(y,j) )^3 - 3*( 1-yTilde(y,j) )*xTilde(x,i)^2 + xTilde(x,i)^3) * zData[i,j] +
        -deltaY[j] * (-( 1-yTilde(y,j) )^2 + ( 1-yTilde(y,j) )^3 + 0.5*( 1-yTilde(y,j) )* xTilde(x,i)^2) * -by[i,j] +
        deltaX[i] * (( 1-yTilde(y,j) ) * xTilde(x,i) - 0.5*xTilde(x,i)^2 - ( 1-yTilde(y,j) )*xTilde(x,i)^2 + 0.5*xTilde(x,i)^3) * bx[i,j] +
        (3*xTilde(x,i)^2 - 3*( 1-yTilde(y,j) ) * xTilde(x,i)^2 - xTilde(x,i)^3)*zData[i+1,j+1] + -deltaY[j]*(0.5*xTilde(x,i)^2-0.5*( 1-yTilde(y,j) )*xTilde(x,i))*-by[i+1,j+1] +
        deltaX[i] * (-xTilde(x,i)^2 + ( 1-yTilde(y,j) )*xTilde(x,i)^2 + 0.5*xTilde(x,i)^3)*bx[i+1,j+1] + (3*( 1-yTilde(y,j) )*xTilde(x,i)^2 - xTilde(x,i)^3)*zData[i+1,j] +
        -deltaY[j] * (-0.5( 1-yTilde(y,j) )*xTilde(x,i)^2)* -by[i+1,j] + deltaX[i]*(-( 1-yTilde(y,j) )*xTilde(x,i)^2 + 0.5*xTilde(x,i)^3)*bx[i+1,j]


    # N = 100 # The sample size in each small square is 4N^2
    @objective(model, Min, sum( sum( 1/(N^2)(
            sum( sum( gamma( 1,x[k],y[l],bx[k,l],by[k,l],i,j )  
                    + gamma( 2,x[k],y[l],bx[k,l],by[k,l],i,j )  
                    + gamma( 3,x[k],y[l],bx[k,l],by[k,l],i,j )  
                    + gamma( 4,x[k],y[l],bx[k,l],by[k,l],i,j )  
                for l in 1:k) for k in 1:N)
        +   sum( sum( gamma( 1,x[k],y[l],bx[k,l],by[k,l],i,j ) 
                    + gamma( 2,x[k],y[l],bx[k,l],by[k,l],i,j ) 
                    + gamma( 3,x[k],y[l],bx[k,l],by[k,l],i,j ) 
                    + gamma( 4,x[k],y[l],bx[k,l],by[k,l],i,j )  
                for l in N:2*N-k) for k in N+1:2*N)
        ) for j in 1:J) for i in 1:I))
    return
end