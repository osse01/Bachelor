include("read_from_txt.jl")

xdata_1 = [0,1,2,3,4,5,6]
ydata_1 = [0,1,2,3,4,5,6]
zdata_1 = [0 0 0 0 0 0 0;
           0 0 0 0 0 0 0;
           0 0 1 1 1 0 0;
           0 0 1 1 1 0 0;
           0 0 1 1 1 0 0;
           0 0 0 0 0 0 0;
           0 0 0 0 0 0 0]
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
# ----------------------------------------------------------------
xdata_4 = [0,1,2,3,4,5,6,7,8,9]
ydata_4 = [0,1,2,3,4,5,6,7,8,9]
zdata_4 = [1 1 1 0 0 4 4 4 4 1;
           1 1 1 0 0 1 1 1 1 1;
           1 1 1 3 3 3 3 3 3 3;
           2 2 2 1 0 0 0 0 0 0;
           3 3 3 1 0 0 0 0 0 0;
           4 4 4 1 0 0 2 2 0 0;
           1 0 3 1 0 0 2 2 0 0;
           1 0 3 1 0 0 0 0 0 0;
           1 0 3 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0 0 0]
# ----------------------------------------------------------------
xdata_5 = [0,1,2,3,4,5,6]
ydata_5 = [0,1,2,3,4,5,6]
zdata_5 = repeat(sin.(xdata_5.*2pi./3),inner=(1,7))
# ----------------------------------------------------------------
data = read_values_from_file("rhine_data_2d_20.txt")
xdata_6 = data.x
ydata_6 = reverse(data.y)
zdata_6 = reverse(data.z, dims=1)
# ----------------------------------------------------------------
data = read_values_from_file("seaside_oregon.txt")
xdata_7 = data.x
ydata_7 = data.y
zdata_7 = data.z