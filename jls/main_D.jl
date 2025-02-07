# find the spectrum of Landau levels in a cosine potential in 2D as a function of flux
# find the corresponding Wannier plot and color gaps in the spectrum according to Chern number

start_time_init = time();

include(joinpath(@__DIR__,"params.jl"))
using .Params                       # interpret parameters from ARGS; phys consts; set up lists
include(joinpath(@__DIR__,"hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis
include(joinpath(@__DIR__,"states.jl"))
using .States                       # wavefunctions and operations for them + density plotting
include(joinpath(@__DIR__,"plotting.jl"))
using .Plt                          # functions for different kinds of plots      

using ProgressMeter
using LinearAlgebra
using Plots

plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local"

args = ARGS
args = ["[101, 150, 0.01, 50, 5, 1]"]

# get parameters from ARGS
p, q, U0, a_in_angstr, NLL, np = Params.parse_arguments_D(args)
a = a_in_angstr * 1e-10                     # lattice const in meters
phi = p/q                                   # unit flux per unit cell
xi0 = sqrt(2π / phi)

# get lists of ky values to iterate over
Nky = 10;                                   # number of ky* points; independent calculations; variation on scale of U0
ky_list = Params.get_ky_list(a, Nky)

# empty lists to store coordinates for plot
energies = Float64[];                       # all energies 
vectors = Vector{ComplexF64}[];                           # repeating phis, for easier plotting later on

# message for time it took to initialise
end_time_init = time();
elapsed_time_init = round(end_time_init - start_time_init; digits = 3);
println("Initialisation done in $elapsed_time_init seconds.")
println("\n\n==================
Calculating the total electronic density of the first $(NLL+1) Landau levels
in a 2D cos potential with strength U=$U0 eV and lattice constant a=$a_in_angstr Å
at flux $phi with particle density $np. Resolution of the calculation is p=$p.\n")


# iterate over lists and diagonalise hamiltonians
start_time_diag = time();                   # set up a clock to monitor elapsed time

@showprogress for ky in ky_list

    H = Hamil.get_full_ham(xi0, ky, U0, a, p, NLL)
    evalsH, evecsH = eigen(H)
    evals_cut, evecs_cut = States.discard_high_energies(evalsH, evecsH, phi, NLL, np)

    energies = [energies; evals_cut]

    for j in eachindex(evals_cut)
        push!(vectors, evecs_cut[:,j])
    end

end

end_time_diag = time();
elapsed_time_diag = round(end_time_diag - start_time_diag; digits = 3);
println("Spectrum has been calculated in $elapsed_time_diag seconds.")



# =============================== PLOTTING ===============================
start_time_plot = time();



# save plots


# end_time_plot = time();
# elapsed_time_plot = round(end_time_plot - start_time_plot; digits = 3);
# println("All plotting done in $elapsed_time_plot seconds.")
# elapsed_time_all = round(end_time_plot - start_time_init; digits = 3);
# println("Code finished running in $elapsed_time_all seconds. Output files can be found in $plot_save_folder_path.")


