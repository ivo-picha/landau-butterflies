# plot a density plot of the electronic density from a pre-saved .npz file

include(joinpath(dirname(@__DIR__), "mods/plotting.jl"))
using .Plt
include(joinpath(dirname(@__DIR__), "mods/params.jl"))
using .Params

using NPZ
using Plots

npz_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/local/dens_grids_p1-q2-U0.1-a50.0-N5-n1.0.npz"
plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local"

data = npzread(npz_path)
data_x = data["x"]
data_y = data["y"]
data_z = data["z"]

npz_name = basename(npz_path)
a_in_angstr = Params.get_a_value(npz_name)
a = a_in_angstr * 1e-10

# generate plot
plot_d = Plt.plot_density(xgrid, ygrid, density_grid, a)
plots_title = string("ϕ=$p/$q, nₚ=$np, U₀=$U0 eV,  a=$a_in_angstr Å,  Nₗₗ=$NLL") # add title to plot
title!(plot_d, plots_title)

# save plot
savefig(plot_d, joinpath(data_save_folder_path, "D_p$p-q$q-U$U0-a$a_in_angstr-N$NLL-n$np.png"))

