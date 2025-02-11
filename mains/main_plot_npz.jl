# plot a density plot of the electronic density from a pre-saved .npz file

include(joinpath(dirname(@__DIR__), "mods/plotting.jl"))
using .Plt

using NPZ
using Plots

npz_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/local/dens_grids_101-100-0.01-50.0-4-1.0.npz"
plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local"

data = npzread(npz_path)
data_x = data["x"]
data_y = data["y"]
data_z = data["z"]

heatmap(data_x, data_y, data_z)


# Plt.plot_density_smooth(data_x, data_y, data_z)

