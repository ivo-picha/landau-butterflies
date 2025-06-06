# find the electronic densities of eigenstates of Landau levels in a cosine potential in 2D as a function of flux


start_time_init = time();

include(joinpath(dirname(@__DIR__),"mods/params.jl"));
using .Params;                       # interpret parameters from ARGS; phys consts; set up lists
include(joinpath(dirname(@__DIR__),"mods/hamiltonian.jl"));
using .Hamil;                        # build a Hamiltonian matrix in a Landau level basis                   # wavefunctions and operations for them + density plotting
include(joinpath(dirname(@__DIR__),"mods/plotting.jl"));
using .Plt;                          # functions for different kinds of plots      
include(joinpath(dirname(@__DIR__),"mods/densities.jl"));
using .Dens;


using ProgressMeter
using LinearAlgebra
using Plots
using NPZ
using Statistics: mean

args = [1, 2, 0.005, 50, 9]

# get parameters from ARGS
p, q, U0, a_in_angstr, NLL = args                     
a = a_in_angstr * 1e-10                     # lattice const in meters
phi = p/q                                   # unit flux per unit cell
xi0 = sqrt(2π / phi)

p=Int(p)
q=Int(q)
NLL=Int(NLL)

# get lists of ky0 and Y values to iterate over
Nky = 17;                                   # number of ky* points; independent calculations; variation on scale of U0
ky_list = Params.get_ky_list(a, Nky)
NY = 17;
Y_list = Params.get_Y_list(NY)

# message for time it took to initialise
end_time_init = time();
elapsed_time_init = round(end_time_init - start_time_init; digits = 3);
println("Initialisation done in $elapsed_time_init seconds.")

println("Calculating spectrum for flux p/q=$(p)/$(q) at U₀=$(U0) and a=$a including the first $(NLL+1) Landau levels.") 


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
sort!(states_vec, by = first)

end_time_diag = time();
elapsed_time_diag = round(end_time_diag - start_time_diag; digits = 3);
println("\nSpectrum has been calculated and cut off in $elapsed_time_diag seconds. Calculating el. densities...")


# =============================== PLOTTING ===============================
start_time_plot = time();

# create a grid and plot the spectrum point by point
N_uc_x = 5                 # number of unit lengths to be plotted in x and y
N_uc_y = N_uc_x
Ngrid = 128                # number of points in each dimension
xplotrange = range(0,N_uc_x*a,Ngrid)
yplotrange = range(0,N_uc_y*a,Ngrid)
xyplotlist = reshape(collect(Iterators.product(xplotrange,yplotrange)),:)

nmlist = [(n,m) for n = 0:NLL for m = 0:(p-1)]

nn=1;

println("Plotting density of state $nn:\n ky0=$(states_vec[nn][2]), Y=$(states_vec[nn][3]), energy $(states_vec[nn][1])eV")

dens_n = Dens.get_density_list(xyplotlist,states_vec[nn],nmlist,phi,a,p);
dens_N = dens_n .* (1/mean(dens_n)); # normalize

# generate plot
plot_d = Plt.plot_density(collect(xplotrange), collect(yplotrange), Float64.(transpose(reshape(dens_N,Ngrid,Ngrid))), a, 1.);
plots_title = string("ϕ=$p/$q, U₀=$U0 eV,  a=$a_in_angstr Å,  Nₗₗ=$NLL, state #$nn"); # add title to plot
title!(plot_d, plots_title)

# save plot
# savefig(plot_d, joinpath(plot_save_folder_path, "DsmT_p$p-q$q-U$U0-a$a_in_angstr-N$NLL-n$np-T$TK.png"))

# end_time_plot = time();
# elapsed_time_plot = round(end_time_plot - start_time_plot; digits = 3);
# println("El. densities calculated in $elapsed_time_plot seconds on a $(Ngrid)x$(Ngrid) grid.")
# elapsed_time_all = round(end_time_plot - start_time_init; digits = 3);
# println("Code finished running in $elapsed_time_all seconds. Output files can be found in $plot_save_folder_path for plots and $data_save_folder_path for npz files.")



# sort(eigvals(Hamil.get_full_ham(sqrt(2π), 0.2*2π/a, 0.5*(2π), 0.2, a, 1, 1)))
# Hamil.matA(0, 2π, 0.0*2π/a, 0.0, 0.1, a, 1)