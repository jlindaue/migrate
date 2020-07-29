using DelimitedFiles
arr=readdlm("trajectory.txt", Float64)
using Gnuplot
# using DataFrames
#plot(select(arr, :column1), select(arr, :column2))
plot(arr[:,1],arr[:,2])
