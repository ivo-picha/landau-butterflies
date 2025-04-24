# plot a density plot of the electronic density from a pre-saved .npz file

include(joinpath(dirname(@__DIR__), "mods/plotting.jl"))
using .Plt
include(joinpath(dirname(@__DIR__), "mods/params.jl"))
using .Params

using NPZ
using Plots

npz_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/local/dens_grids_p1-q2-U0.1-a50.0-N5-n1.0.npz"
plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/npzdens"


function plot_npz(npz_path::String)
    data = npzread(npz_path)
    data_x = data["x"]
    data_y = data["y"]
    data_z = data["z"]

    npz_name = basename(npz_path)
    pd = Params.extract_params(npz_name)

    p = pd[:p]
    q = pd[:q]
    np = parse(Float64, pd[:np])
    U0 = pd[:U0]
    a_in_angstr = pd[:a_in_angstr]
    a = parse(Float64,a_in_angstr) * 1e-10
    NLL = pd[:NLL]
    TK = pd[:TK]

    # generate plot
    plot_d = Plt.plot_density(data_x, data_y, data_z, a, np)
    plots_title = string("ϕ=$p/$q, nₚ=$np, U₀=$U0 eV,  a=$a_in_angstr Å,  Nₗₗ=$NLL, T=$(TK) K") # add title to plot
    title!(plot_d, plots_title)

    # save plot
    savefig(plot_d, joinpath(plot_save_folder_path, "Dnpz_p$p-q$q-U$U0-a$a_in_angstr-N$NLL-n$np-T$TK.png"))    
end

dir = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/mafalda/densities"

for file in readdir(dir)
    plot_npz(joinpath(dir, file))
end


