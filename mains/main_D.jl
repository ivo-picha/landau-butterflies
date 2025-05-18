# find the spectrum of Landau levels in a cosine potential in 2D as a function of flux
# find the corresponding Wannier plot and color gaps in the spectrum according to Chern number

start_time_init = time();

include(joinpath(dirname(@__DIR__),"mods/params.jl"))
using .Params                       # interpret parameters from ARGS; phys consts; set up lists
include(joinpath(dirname(@__DIR__),"mods/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis                   # wavefunctions and operations for them + density plotting
include(joinpath(dirname(@__DIR__),"mods/plotting.jl"))
using .Plt                          # functions for different kinds of plots      
include(joinpath(dirname(@__DIR__),"mods/densities.jl"))
using .Dens


using ProgressMeter
using LinearAlgebra
using Plots
using NPZ
using Statistics: mean

#plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/densities"
#data_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/local/N1/varyphi_U0.01"
plot_save_folder_path = "/users/ivoga/lh/plts/densities"
data_save_folder_path = "/users/ivoga/lh/data"

args = ARGS
#args = ["[1, 2, 0.029, 50, 20, 1., 1]"]

# get parameters from ARGS
p, q, U0, a_in_angstr, NLL, np, TK = Params.parse_arguments_D(args)
a = a_in_angstr * 1e-10                     # lattice const in meters
TeV = TK * Params.kB                        # temperature in eV
phi = p/q                                   # unit flux per unit cell
xi0 = sqrt(2π / phi)

# get lists of ky0 and Y values to iterate over
Nky = 33;                                   # number of ky* points; independent calculations; variation on scale of U0
ky_list = Params.get_ky_list(a, Nky)
NY = 33;
Y_list = Params.get_Y_list(NY)

# message for time it took to initialise
end_time_init = time();
elapsed_time_init = round(end_time_init - start_time_init; digits = 3);
println("Initialisation done in $elapsed_time_init seconds.")

Params.startmessage_Dsm(p, phi, U0, a_in_angstr, NLL, np, TK)


# iterate over lists and diagonalise hamiltonians
start_time_diag = time();                   # set up a clock to monitor elapsed time

# empty list to store information about each state in tuple (eval, ky, evec)
states_vec = Tuple{Float64, Float64, Float64, Vector{ComplexF64}}[];

@showprogress for Y in Y_list
    for ky0 in ky_list

        H = Hamil.get_full_ham(xi0, ky0, Y, U0, a, p, NLL)
        evalsH, evecsH = eigen(H)

        for j in eachindex(evalsH)
            global states_vec
            push!(states_vec, (evalsH[j], ky0, Y, evecsH[:,j]))
        end

    end
end

# sort all states by energy
states_vec_sorted = sort(states_vec, by = first)
states_vec_cut, EF = Dens.discard_high_energies_smear(states_vec_sorted, phi, NLL, np, TeV) # EF is Fermi energy

end_time_diag = time();
elapsed_time_diag = round(end_time_diag - start_time_diag; digits = 3);
println("\nSpectrum has been calculated and cut off in $elapsed_time_diag seconds. Calculating el. densities...")


# =============================== PLOTTING ===============================
start_time_plot = time();

# create a grid and plot the spectrum point by point
N_uc_x = 1                 # number of unit lengths to be plotted in x and y
N_uc_y = N_uc_x
Ngrid = 16                  # number of points in each dimension
xplotrange = range(0,N_uc_x*a,Ngrid)
yplotrange = range(0,N_uc_y*a,Ngrid)
xyplotlist = reshape(collect(Iterators.product(xplotrange,yplotrange)),:)

nmlist = [(n,m) for n = 0:NLL for m = 0:(p-1)]

tot_dens = zeros(Float64,length(xyplotlist))
@showprogress for n in eachindex(states_vec_cut)
    dens_n = Dens.get_density_list(xyplotlist,states_vec_cut[n],nmlist,phi,a,p)
    dens_n_FD = dens_n.*Dens.fermi_dirac(states_vec_cut[n][1],EF,TeV) # add Fermi-Dirac smearing
    global tot_dens
    tot_dens = tot_dens .+ dens_n_FD
end
tot_dens_N = tot_dens .* (np/mean(tot_dens)) # normalize

x_grid = collect(xplotrange)
y_grid = collect(yplotrange)
dens_grid = Float64.(transpose(reshape(tot_dens_N,Ngrid,Ngrid)))

# save data so it can be accessed later; npz format, readable by python as well
npzwrite(joinpath(data_save_folder_path, "dens_grids_p$p-q$q-U$U0-a$a_in_angstr-N$NLL-n$np-T$TK.npz"),
        Dict("x" => x_grid, "y" => y_grid, "z" => dens_grid)
)

# interpolate data for smoother plotting
x_grid_new, y_grid_new, dens_grid_new = Dens.square_interpolate_d(x_grid, y_grid, dens_grid)

# analyze density, get prominence and width of the density in x and y
av_prom_xy, std_prom_xy, av_width_xy, std_width_xy = Dens.get_proms_widths(x_grid_new ./a, y_grid_new ./a, dens_grid_new, N_uc_x)

# generate plot
plot_d = Plt.plot_density(x_grid_new, y_grid_new, dens_grid_new, a, np)
plots_title = string("ϕ=$p/$q, nₚ=$np, U₀=$U0 eV,  a=$a_in_angstr Å,  Nₗₗ=$NLL, T=$TK K \n 
     E=$(round(EF;digits=3)) eV, w/a = $(round(av_width_xy; digits = 3)), pr = $(round(av_prom_xy;digits=3))") # add title to plot
title!(plot_d, plots_title)

# save plot
savefig(plot_d, joinpath(plot_save_folder_path, "D_p$p-q$q-U$U0-a$a_in_angstr-N$NLL-n$np-T$TK.png"))

end_time_plot = time();
elapsed_time_plot = round(end_time_plot - start_time_plot; digits = 3);
println("El. densities calculated in $elapsed_time_plot seconds on a $(Ngrid)x$(Ngrid) grid.")
elapsed_time_all = round(end_time_plot - start_time_init; digits = 3);
println("Code finished running in $elapsed_time_all seconds. Output files can be found in $plot_save_folder_path for plots and $data_save_folder_path for npz files.")
