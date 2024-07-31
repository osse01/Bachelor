
using JuMP, HiGHS
model = Model(HiGHS.Optimizer)

 set_attribute(model, "presolve", "on")
 set_attribute(model, "solver", "ipm") # Interior Point Method

function biCubicSpline(xData, yData, zData, N, lambda)
    I = length(xData)                               # length of xData
    J = length(yData)                               # length of yData
    deltaX = [xData[i+1]-xData[i] for i in 1:I-1]   # x-step length
    deltaY = [yData[j+1]-yData[j] for j in 1:J-1]   # y-step length

    @variable(model, bx[1:I,1:J])               # dzdx(x_i,y_i)
    @variable(model, by[1:I,1:J])               # dzdy(x_i,y_i)
    
    d2zdx2_1(i,j,xTilde, yTilde) = 1/(deltaX[i]^2)*((-6+12*xTilde)*zData[i,j]
                                    +deltaX[i]*(-4+6*xTilde)*bx[i,j]
                                    +(6-12*xTilde)*zData[i+1,j]
                                    +deltaX[i]*(-1+6*deltaX[i])*bx[i+1,j])
    d2z1dxdy_1(i,j,xTilde, yTilde) = 1/(deltaX[i]*deltaY[j])*(
                                    6*yTilde*( (zData[i,j]+zData[i+1,j+1]) - (zData[i+1,j]+zData[i,j+1]) )
                                    +deltaX[i]*yTilde*( (bx[i,j]+bx[i+1,j]) - (bx[i,j+1]+bx[i+1,j+1]) )
                                    +deltaY[j]*( (by[i+1,j]-by[i,j]) + 
                                    2*yTilde*( (by[i,j]+by[i,j+1]) - (by[i+1,j]+by[i+1,j+1]) ) ) )
    d2z1dy2_1(i,j,xTilde,yTilde) = 1/(deltaY[j]^2)*((-6+6*xTilde+6*yTilde)*zData[i,j]
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
                                    +deltaY[j]*(-2*xTilde+3*yTilde)*by[i+1,j+1]
                                    )

    d2zdx2_2(i,j,xTilde, yTilde) = 1/(deltaY[i]^2)*((-6+12*yTilde)*zData[i+1,j]
                                    +deltaY[i]*(-4+6*yTilde)*by[i+1,j]
                                    +(6-12*yTilde)*zData[i+1,j+1]
                                    +deltaY[i]*(-1+6*deltaY[i])*by[i+1,j+1])
    d2z1dxdy_2(i,j,xTilde, yTilde) = 1/(deltaY[i]*-deltaX[j])*(
                                    6*(1-xTilde)*( (zData[i+1,j]+zData[i,j+1]) - (zData[i+1,j+1]+zData[i,j]) )
                                    +deltaY[i]*(1-xTilde)*( (by[i+1,j]+by[i+1,j+1]) - (by[i,j]+by[i,j+1]) )
                                    +(-deltaX[j])*( (-bx[i+1,j+1]+bx[i+1,j]) + 
                                    2*(1-xTilde)*( (bx[i+1,j]+bx[i,j]) - (-bx[i+1,j+1]-bx[i,j+1]) ) ) )
    d2z1dy2_2(i,j,xTilde,yTilde) = 1/(deltaY[i]^2)*((-6+6*yTilde+6*(1-xTilde))*zData[i+1,j]
                                    +deltaY[i]*(-1+yTilde)*by[i+1,j]
                                    -deltaX[j]*(-3+2*yTilde+3*(1-xTilde))*-bx[i+1,j]
                                    +(-6*yTilde+6*(1-xTilde))*zData[i+1,j+1]
                                    +deltaY[i]*(yTilde)*by[i+1,j+1]
                                    -deltaX[j]*(-1-2*yTilde+3*(1-xTilde))*-bx[i+1,j+1]
                                    +(6-6*yTilde-6*(1-xTilde))*zData[i,j]
                                    +deltaY[i]*(1-yTilde)*by[i,j]
                                    -deltaX[j]*(-2+2*yTilde+3*(1-xTilde))*-bx[i,j]
                                    +(6*yTilde-6*(1-xTilde))*zData[i,j+1]
                                    +deltaY[i]*(-yTilde)*by[i,j+1]
                                    -deltaX[j]*(-2*yTilde+3*(1-xTilde))*-bx[i,j+1]
                                    )
    d2zdx2_3(i,j,xTilde, yTilde) = 1/(deltaX[i]^2)*((-6+12*( 1-xTilde ))*zData[i+1,j+1]
                                    -deltaX[i]*(-4+6*( 1-xTilde ))*bx[i+1,j+1]
                                    +(6-12*( 1-xTilde ))*zData[i,j+1]
                                    -deltaX[i]*(-1+6*-deltaX[i])*bx[i,j+1])
    d2z1dxdy_3(i,j,xTilde, yTilde) = 1/(deltaX[i]*deltaY[j])*(
                                    6*( 1-yTilde )*( (zData[i+1,j+1]+zData[i,j]) - (zData[i,j+1]+zData[i+1,j]) )
                                    -deltaX[i]*( 1-yTilde )*( (bx[i+1,j+1]+bx[i,j+1]) - (bx[i+1,j]+bx[i,j]) )
                                    -deltaY[j]*( (by[i,j+1]-by[i+1,j+1]) + 
                                    2*( 1-yTilde )*( (by[i+1,j+1]+by[i+1,j]) - (by[i,j+1]+by[i,j]) ) ) )
    d2z1dy2_3(i,j,xTilde,yTilde) = 1/(deltaY[i]^2)*((-6+6*( 1-xTilde )+6*( 1-yTilde ))*zData[i+1,j+1]
                                    -deltaX[i]*(-1+( 1-xTilde ))*bx[i+1,j+1]
                                    -deltaY[j]*(-3+2*( 1-xTilde )+3*( 1-yTilde ))*by[i+1,j+1]
                                    +(-6*( 1-xTilde )+6*( 1-yTilde ))*zData[i,j+1]
                                    -deltaX[i]*(( 1-xTilde ))*bx[i,j+1]
                                    -deltaY[j]*(-1-2*( 1-xTilde )+3*( 1-yTilde ))*by[i,j+1]
                                    +(6-6*( 1-xTilde )-6*( 1-yTilde ))*zData[i+1,j]
                                    -deltaX[i]*(1-( 1-xTilde ))*bx[i+1,j]
                                    -deltaY[j]*(-2+2*( 1-xTilde )+3*( 1-yTilde ))*by[i+1,j]
                                    +(6*( 1-xTilde )-6*( 1-yTilde ))*zData[i,j]
                                    -deltaX[i]*(-( 1-xTilde ))*bx[i,j]
                                    -deltaY[j]*(-2*( 1-xTilde )+3*( 1-yTilde ))by[i,j]
                                    )
    d2zdx2_4(i,j,xTilde, yTilde) = 1/(deltaY[j]^2)*((-6+12*( 1-yTilde ))*zData[i,j+1]
                                    -deltaY[j]*(-4+6*( 1-yTilde ))*-bx[i,j+1]
                                    +(6-12*( 1-yTilde ))*zData[i,j]
                                    -deltaY[j]*(-1+6*-deltaY[j])*-bx[i,j])
    d2z1dxdy_4(i,j,xTilde, yTilde) = -1/(deltaX[i]*deltaY[j])*(
                                    6*xTilde*( (zData[i,j+1]+zData[i+1,j]) - (zData[i,j]+zData[i+1,j+1]) )
                                    -deltaY[j]*xTilde*( (-bx[i,j+1]+-bx[i,j]) - (-bx[i+1,j+1]+-bx[i+1,j]) )
                                    +deltaX[i]*( (by[i,j]-by[i,j+1]) + 
                                    2*xTilde*( (by[i,j+1]+by[i+1,j+1]) - (by[i,j]+by[i+1,j]) ) ) )
    d2z1dy2_4(i,j,xTilde,yTilde) = 1/(deltaX[i]^2)*((-6+6*( 1-yTilde )+6*xTilde)*zData[i,j+1]
                                    -deltaY[j]*(-1+( 1-yTilde ))*-bx[i,j+1]
                                    +deltaX[i]*(-3+2*( 1-yTilde )+3*xTilde)*by[i,j+1]
                                    +(-6*( 1-yTilde )+6*xTilde)*zData[i,j]
                                    -deltaY[j]*(( 1-yTilde ))*-bx[i,j]
                                    +deltaX[i]*(-1-2*( 1-yTilde )+3*xTilde)*by[i,j]
                                    +(6-6*( 1-yTilde )-6*xTilde)*zData[i+1,j+1]
                                    -deltaY[j]*(1-( 1-yTilde ))*-bx[i+1,j+1]
                                    +deltaX[i]*(-2+2*( 1-yTilde )+3*yTilde)*by[i+1,j+1]
                                    +(6*( 1-yTilde )-6*xTilde)*zData[i+1,j]
                                    -deltaY[j]*(-( 1-yTilde ))*-bx[i+1,j]
                                    +deltaX[i]*(-2*( 1-yTilde )+3*xTilde)by[i+1,j]
                                    )

    @variable(model, abs_d2zdx2_1[1:I-1,1:J-1, 1:(2*N), 1:(2*N) ]>=0)
    @variable(model, abs_d2z1dxdy_1[1:I-1,1:J-1, 1:(2*N), 1:(2*N) ]>=0)
    @variable(model, abs_d2z1dy2_1[1:I-1,1:J-1, 1:(2*N), 1:(2*N) ]>=0)
    @variable(model, abs_d2zdx2_2[1:I-1,1:J-1, 1:(2*N), 1:(2*N) ]>=0)
    @variable(model, abs_d2z1dxdy_2[1:I-1,1:J-1, 1:(2*N), 1:(2*N) ]>=0)
    @variable(model, abs_d2z1dy2_2[1:I-1,1:J-1, 1:(2*N), 1:(2*N) ]>=0)
    @variable(model, abs_d2zdx2_3[1:I-1,1:J-1, 1:(2*N), 1:(2*N) ]>=0)
    @variable(model, abs_d2z1dxdy_3[1:I-1,1:J-1, 1:(2*N), 1:(2*N) ]>=0)
    @variable(model, abs_d2z1dy2_3[1:I-1,1:J-1, 1:(2*N), 1:(2*N) ]>=0)
    @variable(model, abs_d2zdx2_4[1:I-1,1:J-1, 1:(2*N), 1:(2*N) ]>=0)
    @variable(model, abs_d2z1dxdy_4[1:I-1,1:J-1, 1:(2*N), 1:(2*N) ]>=0)
    @variable(model, abs_d2z1dy2_4[1:I-1,1:J-1, 1:(2*N), 1:(2*N) ]>=0)
    @variable(model, abs_bx[1:I,1:J]>=0)
    @variable(model, abs_by[1:I,1:J]>=0)

    gamma_1(i, j, k, l ) = 
        abs_d2zdx2_1[i,j, k, l ] + 2*abs_d2z1dxdy_1[i,j, k, l ] + abs_d2z1dy2_1[i,j, k, l ]
    gamma_2(i, j, k, l ) = 
        abs_d2zdx2_2[i,j, k, l ] + 2*abs_d2z1dxdy_2[i,j, k, l ] + abs_d2z1dy2_2[i,j, k, l ]
    gamma_3(i, j, k, l ) = 
        abs_d2zdx2_3[i,j, k, l ] + 2*abs_d2z1dxdy_3[i,j, k, l ] + abs_d2z1dy2_3[i,j, k, l ]
    gamma_4(i, j, k, l ) = 
        abs_d2zdx2_4[i,j, k, l ] + 2*abs_d2z1dxdy_4[i,j, k, l ] + abs_d2z1dy2_4[i,j, k, l ]
    
        # The sample size in each small square is 4N^2
    # lambda is a small number
    @objective(model, Min, sum( sum( 1/(N^2)*(
            sum( sum( abs_d2zdx2_1[i,j, k, l ] + 2*abs_d2z1dxdy_1[i,j, k, l ] + abs_d2z1dy2_1[i,j, k, l ]  
                    + abs_d2zdx2_2[i,j, k, l ] + 2*abs_d2z1dxdy_2[i,j, k, l ] + abs_d2z1dy2_2[i,j, k, l ] 
                    + abs_d2zdx2_3[i,j, k, l ] + 2*abs_d2z1dxdy_3[i,j, k, l ] + abs_d2z1dy2_3[i,j, k, l ] 
                    + abs_d2zdx2_4[i,j, k, l ] + 2*abs_d2z1dxdy_4[i,j, k, l ] + abs_d2z1dy2_4[i,j, k, l ] 
                for l in 1:k) for k in 1:N)
        +   sum( sum( abs_d2zdx2_1[i,j, k, l ] + 2*abs_d2z1dxdy_1[i,j, k, l ] + abs_d2z1dy2_1[i,j, k, l ]  
                    + abs_d2zdx2_2[i,j, k, l ] + 2*abs_d2z1dxdy_2[i,j, k, l ] + abs_d2z1dy2_2[i,j, k, l ] 
                    + abs_d2zdx2_3[i,j, k, l ] + 2*abs_d2z1dxdy_3[i,j, k, l ] + abs_d2z1dy2_3[i,j, k, l ] 
                    + abs_d2zdx2_4[i,j, k, l ] + 2*abs_d2z1dxdy_4[i,j, k, l ] + abs_d2z1dy2_4[i,j, k, l ]   
                for l in N:2*N-k) for k in N+1:2*N)
        ) for j in 1:J-1) for i in 1:I-1)) + 
            lambda*sum(sum( abs_bx[i,j] + abs_by[i,j] for i in 1:I) for j in 1:J)

    @constraint(model,pos_d2zdx2_1[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2zdx2_1[i, j, k, l] >= d2zdx2_1(i, j, k/(2*N), l/(2*N)) )
    @constraint(model,neg_d2zdx2_1[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2zdx2_1[i, j, k, l] >= -d2zdx2_1(i, j, k/(2*N), l/(2*N)))

    @constraint(model,pos_d2z1dxdy_1[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dxdy_1[i, j, k, l] >= d2z1dxdy_1(i, j, k/(2*N), l/(2*N)) )
    @constraint(model,neg_d2z1dxdy_1[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dxdy_1[i, j, k, l] >= -d2z1dxdy_1(i, j, k/(2*N), l/(2*N)))

    @constraint(model,pos_d2z1dy2_1[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dy2_1[i, j, k, l] >= d2z1dy2_1(i, j, k/(2*N), l/(2*N)) )
    @constraint(model,neg_d2z1dy2_1[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dy2_1[i, j, k, l] >= -d2z1dy2_1(i, j, k/(2*N), l/(2*N)))

    @constraint(model,pos_d2zdx2_2[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2zdx2_2[i, j, k, l] >= d2zdx2_2(i, j, k/(2*N), l/(2*N)) )
    @constraint(model,neg_d2zdx2_2[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2zdx2_2[i, j, k, l] >= -d2zdx2_2(i, j, k/(2*N), l/(2*N)))

    @constraint(model,pos_d2z1dxdy_2[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dxdy_2[i, j, k, l] >= d2z1dxdy_2(i, j, k/(2*N), l/(2*N)) )
    @constraint(model,neg_d2z1dxdy_2[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dxdy_2[i, j, k, l] >= -d2z1dxdy_2(i, j, k/(2*N), l/(2*N)))

    @constraint(model,pos_d2z1dy2_2[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dy2_2[i, j, k, l] >= d2z1dy2_2(i, j, k/(2*N), l/(2*N)) )
    @constraint(model,neg_d2z1dy2_2[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dy2_2[i, j, k, l] >= -d2z1dy2_2(i, j, k/(2*N), l/(2*N)))

    @constraint(model,pos_d2zdx2_3[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2zdx2_3[i, j, k, l] >= d2zdx2_3(i, j, k/(2*N), l/(2*N)) )
    @constraint(model,neg_d2zdx2_3[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2zdx2_3[i, j, k, l] >= -d2zdx2_3(i, j, k/(2*N), l/(2*N)))

    @constraint(model,pos_d2z1dxdy_3[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dxdy_3[i, j, k, l] >= d2z1dxdy_3(i, j, k/(2*N), l/(2*N)) )
    @constraint(model,neg_d2z1dxdy_3[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dxdy_3[i, j, k, l] >= -d2z1dxdy_3(i, j, k/(2*N), l/(2*N)))

    @constraint(model,pos_d2z1dy2_3[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dy2_3[i, j, k, l] >= d2z1dy2_3(i, j, k/(2*N), l/(2*N)) )
    @constraint(model,neg_d2z1dy2_3[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dy2_3[i, j, k, l] >= -d2z1dy2_3(i, j, k/(2*N), l/(2*N)))

    @constraint(model,pos_d2zdx2_4[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2zdx2_4[i, j, k, l] >= d2zdx2_4(i, j, k/(2*N), l/(2*N)) )
    @constraint(model,neg_d2zdx2_4[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2zdx2_4[i, j, k, l] >= -d2zdx2_4(i, j, k/(2*N), l/(2*N)))

    @constraint(model,pos_d2z1dxdy_4[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dxdy_4[i, j, k, l] >= d2z1dxdy_4(i, j, k/(2*N), l/(2*N)) )
    @constraint(model,neg_d2z1dxdy_4[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dxdy_4[i, j, k, l] >= -d2z1dxdy_4(i, j, k/(2*N), l/(2*N)))

    @constraint(model,pos_d2z1dy2_4[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dy2_4[i, j, k, l] >= d2z1dy2_4(i, j, k/(2*N), l/(2*N)) )
    @constraint(model,neg_d2z1dy2_4[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N) ], abs_d2z1dy2_4[i, j, k, l] >= -d2z1dy2_4(i, j, k/(2*N), l/(2*N)))

    @constraint(model,pos_bx[i in 1:I-1, j in 1:J-1], abs_bx[i,j] >= bx[i,j])
    @constraint(model,neg_bx[i in 1:I-1, j in 1:J-1], abs_bx[i,j] >= -bx[i,j])

    @constraint(model,pos_by[i in 1:I-1, j in 1:J], abs_by[i,j] >= by[i,j])
    @constraint(model,neg_by[i in 1:I-1, j in 1:J], abs_by[i,j] >= -by[i,j])

    optimize!(model)
    
    bx = Array(value.(bx))
    by = Array(value.(by))

    print(bx)
    print(by)
    print("\n")
    return [xdata, ydata, zdata, bx, by]
end

function evaluate(spline, N, M)
    xData = spline[1]
    yData = spline[2]
    zData = spline[3]
    bx = spline[4]
    by = spline[5]
    I = length(xData)                               # length of xData
    J = length(yData)                               # length of yData
    deltaX = [xData[i+1]-xData[i] for i in 1:I-1]   # x-step length
    deltaY = [yData[j+1]-yData[j] for j in 1:J-1]   # y-step length

    xTilde(x,i) = (x - xData[i]) / deltaX[i]        # Transforms [x_i,x_i+1] to [0,1]
    yTilde(y,j) = (y - yData[j]) / deltaY[j]        # Transforms [y_i,y_i+1] to [0,1]

    z1(x,y,bx,by,i,j) =
        (1 - 3*xTilde(x,i)^2 + 2*xTilde(x,i)^3 - 3*yTilde(y,j)^2 + 3*xTilde(x,i)*yTilde(y,j)^2 + yTilde(y,j)^3)*zData[i,j] +
        deltaX[i] * (xTilde(x,i) - 2*xTilde(x,i)^2 + xTilde(x,i)^3 - 0.5*yTilde(y,j)^2 + 0.5 * xTilde(x,i)*yTilde(y,j)^2) * bx[i,j] + 
        deltaY[j] * (yTilde(y,j) - xTilde(x,i)*yTilde(y,j) - 1.5 * yTilde(y,j)^2 + xTilde(x,i)*yTilde(y,j)^2 + 0.5*yTilde(y,j)^3) * by[i,j] +
        (3*xTilde(x,i)^2 - 2*xTilde(x,i)^3 - 3*xTilde(x,i)*yTilde(y,j)^2 + yTilde(y,j)^3) * zData[i+1,j] +
        deltaX[i] * (-xTilde(x,i)^2 + xTilde(x,i)^3 + 0.5*xTilde(x,i)* yTilde(y,j)^2) * bx[i+1,j] +
        deltaY[j] * (xTilde(x,i) * yTilde(y,j) - 0.5*yTilde(y,j)^2 - xTilde(x,i)*yTilde(y,j)^2 + 0.5*yTilde(y,j)^3) * by[i+1,j] +
        (3*yTilde(y,j)^2 - 3*xTilde(x,i) * yTilde(y,j)^2 - yTilde(y,j)^3)*zData[i,j+1] + deltaX[i]*(0.5*yTilde(y,j)^2-0.5*xTilde(x,i)*yTilde(y,j)^2)*bx[i,j+1] +
        deltaY[j] * (-yTilde(y,j)^2 + xTilde(x,i)*yTilde(y,j)^2 + 0.5*yTilde(y,j)^3)*by[i,j+1] + (3*xTilde(x,i)*yTilde(y,j)^2 - yTilde(y,j)^3)*zData[i+1,j+1] +
        deltaX[i] * (-0.5*xTilde(x,i)*yTilde(y,j)^2)* bx[i+1,j+1] + deltaY[j]*(-xTilde(x,i)*yTilde(y,j)^2 + 0.5*yTilde(y,j)^3)*by[i+1,j+1]

    z2(x,y,bx,by,i,j) =
        (1 - 3*yTilde(y,j)^2 + 2*yTilde(y,j)^3 - 3*( 1-xTilde(x,i) )^2 + 3*yTilde(y,j)*( 1-xTilde(x,i) )^2 + ( 1-xTilde(x,i) )^3)*zData[i+1,j] +
        deltaY[j] * (yTilde(y,j) - 2*yTilde(y,j)^2 + yTilde(y,j)^3 - 0.5*(1-xTilde(x,i))^2 + 0.5 * yTilde(y,j)*( 1-xTilde(x,i) )^2) * by[i+1,j] + 
        -deltaX[i] * (( 1-xTilde(x,i) ) - yTilde(y,j)*( 1-xTilde(x,i) ) - 1.5 * ( 1-xTilde(x,i) )^2 + yTilde(y,j)*( 1-xTilde(x,i) )^2 + 0.5*( 1-xTilde(x,i) )^3) * -bx[i+1,j] +
        (3*yTilde(y,j)^2 - 2*yTilde(y,j)^3 - 3*yTilde(y,j)*( 1-xTilde(x,i) )^2 + ( 1-xTilde(x,i) )^3) * zData[i+1,j+1] +
        deltaY[j] * (-yTilde(y,j)^2 + yTilde(y,j)^3 + 0.5*yTilde(y,j)* ( 1-xTilde(x,i) )^2) * by[i+1,j+1] +
        -deltaX[i] * (yTilde(y,j) * ( 1-xTilde(x,i) ) - 0.5*( 1-xTilde(x,i) )^2 - yTilde(y,j)*( 1-xTilde(x,i) )^2 + 0.5*( 1-xTilde(x,i) )^3) * -bx[i+1,j+1] +
        (3*( 1-xTilde(x,i) )^2 - 3*yTilde(y,j) * ( 1-xTilde(x,i) )^2 - ( 1-xTilde(x,i) )^3)*zData[i,j] + deltaY[j]*(0.5*( 1-xTilde(x,i) )^2-0.5*yTilde(y,j)*( 1-xTilde(x,i) )^2)*by[i,j] +
        -deltaX[i] * (-( 1-xTilde(x,i) )^2 + yTilde(y,j)*( 1-xTilde(x,i) )^2 + 0.5*( 1-xTilde(x,i) )^3)*-bx[i,j] + (3*yTilde(y,j)*( 1-xTilde(x,i) )^2 - ( 1-xTilde(x,i) )^3)*zData[i,j+1] +
        deltaY[j] * (-0.5*yTilde(y,j)*( 1-xTilde(x,i) )^2)* by[i,j+1] - deltaX[i]*(-yTilde(y,j)*( 1-xTilde(x,i) )^2 + 0.5*( 1-xTilde(x,i) )^3)*-bx[i,j+1]
    
    z3(x,y,bx,by,i,j) =
        (1 - 3*( 1-xTilde(x,i) )^2 + 2*( 1-xTilde(x,i) )^3 - 3*( 1-yTilde(y,j) )^2 + 3*( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2 + ( 1-yTilde(y,j) )^3)*zData[i+1,j+1] +
        -deltaX[i] * (( 1-xTilde(x,i) ) - 2*( 1-xTilde(x,i) )^2 + ( 1-xTilde(x,i) )^3 - 0.5*( 1-yTilde(y,j) )^2 + 0.5 * ( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2) * -bx[i+1,j+1] + 
        -deltaY[j] * (( 1-yTilde(y,j) ) - ( 1-xTilde(x,i) )*( 1-yTilde(y,j) ) - 1.5 * ( 1-yTilde(y,j) )^2 + ( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2 + 0.5*( 1-yTilde(y,j) )^3) * -by[i+1,j+1] +
        (3*( 1-xTilde(x,i) )^2 - 2*( 1-xTilde(x,i) )^3 - 3*( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2 + ( 1-yTilde(y,j) )^3) * zData[i,j+1] +
        -deltaX[i] * (-( 1-xTilde(x,i) )^2 + ( 1-xTilde(x,i) )^3 + 0.5*( 1-xTilde(x,i) )* ( 1-yTilde(y,j) )^2) * -bx[i,j+1] +
        -deltaY[j] * (( 1-xTilde(x,i) ) * ( 1-yTilde(y,j) ) - 0.5*( 1-yTilde(y,j) )^2 - ( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2 + 0.5*( 1-yTilde(y,j) )^3) * -by[i,j+1] +
        (3*( 1-yTilde(y,j) )^2 - 3*( 1-xTilde(x,i) ) * ( 1-yTilde(y,j) )^2 - ( 1-yTilde(y,j) )^3)*zData[i+1,j] + -deltaX[i]*(0.5*( 1-yTilde(y,j) )^2-0.5*( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2)*-bx[i+1,j] +
        -deltaY[j] * (-( 1-yTilde(y,j) )^2 + ( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2 + 0.5*( 1-yTilde(y,j) )^3)*-by[i+1,j] + (3*( 1-xTilde(x,i) )*( -yTilde(y,j) )^2 - ( 1-yTilde(y,j) )^3)*zData[i,j] +
        -deltaX[i] * (-0.5*( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2)* -bx[i,j] + -deltaY[j]*(-( 1-xTilde(x,i) )*( 1-yTilde(y,j) )^2 + 0.5*( 1-yTilde(y,j) )^3)*-by[i,j]

    z4(x,y,bx,by,i,j) =
        (1 - 3*( 1-yTilde(y,j) )^2 + 2*( 1-yTilde(y,j) )^3 - 3*xTilde(x,i)^2 + 3*( 1-yTilde(y,j) )*xTilde(x,i)^2 + xTilde(x,i)^3)*zData[i,j+1] +
        -deltaY[j] * (( 1-yTilde(y,j) ) - 2*( 1-yTilde(y,j) )^2 + ( 1-yTilde(y,j) )^3 - 0.5*xTilde(x,i)^2 + 0.5 * ( 1-yTilde(y,j) )*xTilde(x,i)^2) * -by[i,j+1] + 
        deltaX[i] * (xTilde(x,i) - ( 1-yTilde(y,j) )*xTilde(x,i) - 1.5 * xTilde(x,i)^2 + ( 1-yTilde(y,j) )*xTilde(x,i)^2 + 0.5*xTilde(x,i)^3) * bx[i,j+1] +
        (3*( 1-yTilde(y,j) )^2 - 2*( 1-yTilde(y,j) )^3 - 3*( 1-yTilde(y,j) )*xTilde(x,i)^2 + xTilde(x,i)^3) * zData[i,j] +
        -deltaY[j] * (-( 1-yTilde(y,j) )^2 + ( 1-yTilde(y,j) )^3 + 0.5*( 1-yTilde(y,j) )* xTilde(x,i)^2) * -by[i,j] +
        deltaX[i] * (( 1-yTilde(y,j) ) * xTilde(x,i) - 0.5*xTilde(x,i)^2 - ( 1-yTilde(y,j) )*xTilde(x,i)^2 + 0.5*xTilde(x,i)^3) * bx[i,j] +
        (3*xTilde(x,i)^2 - 3*( 1-yTilde(y,j) ) * xTilde(x,i)^2 - xTilde(x,i)^3)*zData[i+1,j+1] + -deltaY[j]*(0.5*xTilde(x,i)^2-0.5*( 1-yTilde(y,j) )*xTilde(x,i)^2)*-by[i+1,j+1] +
        deltaX[i] * (-xTilde(x,i)^2 + ( 1-yTilde(y,j) )*xTilde(x,i)^2 + 0.5*xTilde(x,i)^3)*bx[i+1,j+1] + (3*( 1-yTilde(y,j) )*xTilde(x,i)^2 - xTilde(x,i)^3)*zData[i+1,j] +
        -deltaY[j] * (-0.5*( 1-yTilde(y,j) )*xTilde(x,i)^2)* -by[i+1,j] + deltaX[i]*(-( 1-yTilde(y,j) )*xTilde(x,i)^2 + 0.5*xTilde(x,i)^3)*bx[i+1,j]

    #gradient = Array([[deltaY[j]/deltaX[i] for i in 1:I-1] for j in 1:J-1])
    gradient = zeros(I, J)
    for i in 1:I-1
        for j in 1:J-1
            gradient[i,j] = deltaY[j]/deltaX[i]
        end
    end
    
    z = zeros(N,M)
    i = 1
    k = 1
    for x in range(xData[1],xData[end],N)
        j = 1
        l = 1
        for y in range(yData[1],yData[end],M)
            if y <= yData[j] + (x - xData[i]) * gradient[i,j] && y <= yData[j+1] - (x - xData[i]) * gradient[i,j]
                z[k,l] = z1(x,y,bx,by,i,j)
            elseif y >= yData[j] + (x - xData[i]) * gradient[i,j] && y <= yData[j+1] - (x - xData[i]) * gradient[i,j]
                z[k,l] = z4(x,y,bx,by,i,j)
            elseif y >= yData[j+1] - (x - xData[i]) * gradient[i,j] && y >= yData[j] + (x - xData[i]) * gradient[i,j]
                z[k,l] = z3(x,y,bx,by,i,j)
            elseif y <= yData[j] + (x - xData[i]) * gradient[i,j] && y >= yData[j+1] - (x - xData[i]) * gradient[i,j]
                z[k,l] = z2(x,y,bx,by,i,j)
            end
            l = l+1
            if y > yData[j+1]
                j = j+1
            end
        end
        k = k+1
        if x > xData[i+1]
            i = i+1
        end
    end
    return z
end