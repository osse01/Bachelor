

xdata_1 = [0,1,2,3,4,5,6]
ydata_1 = [0,1,2,3,4,5,6]
zdata_1 = [0 0 0 0 0 0 0;
           0 0 0 0 0 0 0;
           0 0 1 1 1 0 0;
           0 0 1 1 1 0 0;
           0 0 1 1 1 0 0;
           0 0 0 0 0 0 0;
           0 0 0 0 0 0 0]

xPlotData_1 = repeat(xdata_1,length(xdata_1))
yPlotData_1 = repeat(xdata_1,inner=length(xdata_1))
zPlotData_1 = vec(zdata_1)
# ----------------------------------------------------------------
xdata_2 = [0,1,2,3,4,5,6]
ydata_2 = [0,1,2,3,4,5,6]
zdata_2 = [0 0 0 1 1 1 1;
           0 0 0 1 1 1 1;
           0 0 0 1 1 1 1;
           0 0 0 1 1 1 1;
           0 0 0 1 1 1 1;
           0 0 0 1 1 1 1;
           0 0 0 1 1 1 1]

xPlotData_2 = repeat(xdata_2,length(xdata_2))
yPlotData_2 = repeat(xdata_2,inner=length(xdata_2))
zPlotData_2 = vec(zdata_2)
# ----------------------------------------------------------------
xdata_3 = [0,1,2,3,4,5,6]
ydata_3 = [0,1,2,3,4,5,6]
zdata_3 = [1 1 1 1 1 1 1;
           1 1 1 1 1 1 1;
           1 1 1 1 1 1 1;
           1 1 1 1 1 1 1;
           1 1 1 1 1 1 1;
           1 1 1 1 1 1 1;
           1 1 1 1 1 1 1]

xPlotData_3 = repeat(xdata_3,length(xdata_3))
yPlotData_3 = repeat(xdata_3,inner=length(xdata_3))
zPlotData_3 = vec(zdata_3)