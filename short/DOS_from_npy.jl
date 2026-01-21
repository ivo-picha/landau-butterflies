using Plots
using NPZ
using KernelDensity

#input folder with npy files
infolder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/data_DOS/"
#output folder
outfolder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/plots_DOS/"
mkpath(outfolder)

plt = plot();
for file in readdir(infolder)
    if endswith(file, ".npy")
        filepath = joinpath(infolder, file)
        data = npzread(filepath)
        eta = (maximum(data)-minimum(data))/100
        kd = kde(data, bandwidth = eta)
        plot!(plt, kd.x, kd.density, linewidth = 2, label = "")
    end
end


plt