using DelimitedFiles
arr=readdlm("trajectory.txt", Float64)
using Plots
display(scatter(arr[1,:],arr[2,:]))

