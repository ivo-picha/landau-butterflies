# find the spectrum of Landau levels in a cosine potential in 2D as a function of flux
# find the corresponding Wannier plot and color gaps in the spectrum according to Chern number

start_time_init = time();

include(joinpath(dirname(@__DIR__),"mods/params.jl"))
using .Params                       # interpret parameters from ARGS; phys consts; set up lists
include(joinpath(dirname(@__DIR__),"mods/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis
include(joinpath(dirname(@__DIR__),"mods/states.jl"))
using .States                       # wavefunctions and operations for them + density plotting
include(joinpath(dirname(@__DIR__),"mods/plotting.jl"))
using .Plt                          # functions for different kinds of plots      

using ProgressMeter
using LinearAlgebra
using Plots
using NPZ

plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local"
data_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/local"

args = ARGS
args = ["[15, 7, -0.1, 50, 10, 1]"]

# get parameters from ARGS
p, q, U0, a_in_angstr, NLL, np = Params.parse_arguments_D(args)
a = a_in_angstr * 1e-10                     # lattice const in meters
phi = p/q                                   # unit flux per unit cell
xi0 = sqrt(2π / phi)

# get lists of ky values to iterate over
Nky = 20;                                   # number of ky* points; independent calculations; variation on scale of U0
ky_list = Params.get_ky_list(a, Nky)

# empty lists to store coordinates for plot
energies = Float64[];                       # all energies 
vectors = Vector{ComplexF64}[];                           # repeating phis, for easier plotting later on

# message for time it took to initialise
end_time_init = time();
elapsed_time_init = round(end_time_init - start_time_init; digits = 3);
println("Initialisation done in $elapsed_time_init seconds.")

Params.startmessage_D(p, phi, U0, a_in_angstr, NLL, np)


# iterate over lists and diagonalise hamiltonians
start_time_diag = time();                   # set up a clock to monitor elapsed time

@showprogress for ky in ky_list

    H = Hamil.get_full_ham(xi0, ky, U0, a, p, NLL)
    evalsH, evecsH = eigen(H)
    evals_cut, evecs_cut = States.discard_high_energies(evalsH, evecsH, phi, NLL, np)

    global energies
    energies = [energies; evals_cut]

    for j in eachindex(evals_cut)
        global vectors
        push!(vectors, evecs_cut[:,j])
    end

end

end_time_diag = time();
elapsed_time_diag = round(end_time_diag - start_time_diag; digits = 3);
println("\nSpectrum has been calculated in $elapsed_time_diag seconds. Calculating el. densities...")



# =============================== PLOTTING ===============================
start_time_plot = time();

# create a grid and plot the spectrum point by point
N_uc_x = 3                 # number of unit lengths to be plotted in x
N_uc_y = 3
Nppuc = 20                  # number of points per unit length
xgrid, ygrid, density_grid = States.get_density_grids(N_uc_x, N_uc_y, Nppuc, vectors, collect(ky_list), phi, a, p, NLL)

# save data so it can be accessed later; npz format, readable by python as well
npzwrite(joinpath(data_save_folder_path, "dens_grids_p$p-q$q-U$U0-a$a_in_angstr-N$NLL-n$np.npz"),
        Dict("x" => collect(xgrid), "y" => collect(ygrid), "z" => Float64.(density_grid))
)

# generate plot
plot_d = Plt.plot_density(xgrid, ygrid, density_grid, a)
plots_title = string("ϕ=$p/$q, nₚ=$np, U₀=$U0 eV,  a=$a_in_angstr Å,  Nₗₗ=$NLL") # add title to plot
title!(plot_d, plots_title)

# save plot
savefig(plot_d, joinpath(plot_save_folder_path, "D_p$p-q$q-U$U0-a$a_in_angstr-N$NLL-n$np.png"))

end_time_plot = time();
elapsed_time_plot = round(end_time_plot - start_time_plot; digits = 3);
println("El. densities calculated in $elapsed_time_plot seconds on a $(N_uc_x*Nppuc)x$(N_uc_y*Nppuc) grid.")
elapsed_time_all = round(end_time_plot - start_time_init; digits = 3);
println("Code finished running in $elapsed_time_all seconds. Output files can be found in $plot_save_folder_path for plots and $data_save_folder_path for npz files.")
