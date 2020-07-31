using DelimitedFiles
arr=readdlm("trajectory.txt", Float64)
using Plots
display(scatter(arr[1,:],arr[2,:]))
#display(scatter(arr[6,:],arr[3,:]))
#display(scatter(arr[6,:],arr[4,:]))
#display(scatter(arr[6,:],arr[5,:]))

