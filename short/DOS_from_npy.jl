using Plots
using NPZ
using KernelDensity

#input folder with npy files
infolder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/data_DOS/"
#output folder
outfolder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/plots_DOS/"
mkpath(outfolder)

kd_list = []
for file in readdir(infolder)
    if endswith(file, ".npy")
        filepath = joinpath(infolder, file)
        data = npzread(filepath)
        eta = (maximum(data)-minimum(data))/100
        kd = kde(data, bandwidth = eta)
        push!(kd_list, kd)
    end
end

readdir(infolder) 

kd_list[1].x
# plotting
plt = Plots.plot(framestyle = :box,
                xaxis = :flip,
                xticks = [],
                yticks = false,
                legend = false,
                size = (1000,200),
                dpi = 900,
                xlims = (-0.0141, -0.0018)

)
lww = 5
plot!(plt, kd_list[end].x, kd_list[end].density/maximum(kd_list[end].density), label="ϕ=0", lw=lww, color = :black)
plot!(plt, kd_list[1].x, kd_list[1].density/maximum(kd_list[1].density), label="ϕ=1", lw=lww, color = :darkorchid2)
plot!(plt, kd_list[2].x, kd_list[2].density/maximum(kd_list[2].density), label="ϕ=1", lw=lww, color = :crimson)
plot!(plt, kd_list[3].x, kd_list[3].density/maximum(kd_list[3].density), label="ϕ=2", lw=lww, color = :darkorange2)

savefig(plt, joinpath(outfolder, "DOS_comparison_mainfig1.png"))
#0.027-(17.8-11.8)/(26.5-17.8)*0.001
#0.033+(76.485-69.6)/(26.5-17.8)*0.001
