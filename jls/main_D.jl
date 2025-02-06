# find the spectrum of Landau levels in a cosine potential in 2D as a function of flux
# find the corresponding Wannier plot and color gaps in the spectrum according to Chern number

start_time_init = time();

include(joinpath(@__DIR__,"params.jl"))
using .Params                       # interpret parameters from ARGS; phys consts; set up lists
include(joinpath(@__DIR__,"hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis

include(joinpath(@__DIR__,"plotting.jl"))
using .Plt                          # functions for different kinds of plots      

using ProgressMeter
using LinearAlgebra
using Plots

plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local"

args = ARGS
args = ["[0.05, 1, 0.01 , 50, 101, 1, 0.4]"]

# get parameters from ARGS
startphi, endphi, U0, a_in_angstr, p, NLL, gap_factor = Params.parse_arguments(args)
a = a_in_angstr * 1e-10                     # lattice const in meters

# get lists of q and ky values to iterate over
q_list = Params.get_q_list(startphi, endphi, p)
Nky = 6;                                    # number of ky* points; independent calculations; variation on scale of U0
ky_list = Params.get_ky_list(a, Nky)

# empty lists to store coordinates for plot
energies = Float64[];                       # all energies 
phis = Float64[];                           # repeating phis, for easier plotting later on
particle_densities_global = Float64[];      # particle densities beneath gaps
phis_w = Float64[];                         # separate list to store coordinates for wannier plots
gaps_global = Float64[];                    # gap sizes for size of Wannier plot points
energies_lower_gaps_global = Float64[];     # indexing at what energy at a given flux the gap appears

# set up a threshold for what is identified as a gap; factor of difference in LLs at start
gap_threshold = gap_factor*(Hamil.E_LL(1,sqrt(2π/startphi),a) - Hamil.E_LL(0,sqrt(2π/startphi),a))

# message for time it took to initialise
end_time_init = time();
elapsed_time_init = round(end_time_init - start_time_init; digits = 3);
println("Initialisation done in $elapsed_time_init seconds.")
println("\n\n==================
Calculating the spectrum of the first $(NLL+1) Landau levels
in a 2D cos potential with strength U=$U0 eV and lattice constant a=$a_in_angstr A
from flux $startphi to flux $endphi. Resolution of the calculation is p=$p.\n")
# message for size of calculation
Params.print_size_message(q_list, p, Nky, NLL)

# iterate over lists and diagonalise hamiltonians
start_time_diag = time();                   # set up a clock to monitor elapsed time
@showprogress for q in q_list

    phi = p/q                               # unit flux per unit cell
    xi0 = sqrt(2π / phi)
    energies_at_phi = Float64[];            # to be appended to global energies list

    for ky in ky_list

        H = Hamil.get_full_ham(xi0, ky, U0, a, p, NLL)
        evalsH = eigvals(H)
        energies_at_phi = [energies_at_phi; evalsH] #add eigenvalues to list of energies

    end

    global energies
    energies = [energies; energies_at_phi]
    global phis
    phis = [phis; [phi for j = 1:length(energies_at_phi)]] #add phi values to list of phis

    # get wannier information
    densities, gaps, energies_lower_gaps = Wannier.get_densities_gaps(energies_at_phi, gap_threshold, phi, NLL)

    global particle_densities_global
    particle_densities_global = [particle_densities_global; densities]
    global phis_w
    phis_w = [phis_w; [phi for j = 1:length(densities)]]
    global gaps_global
    gaps_global = [gaps_global; gaps]
    global energies_lower_gaps_global
    energies_lower_gaps_global = [energies_lower_gaps_global; energies_lower_gaps]

end
end_time_diag = time();
elapsed_time_diag = round(end_time_diag - start_time_diag; digits = 3);
println("Spectrum and Wannier plot have been calculated in $elapsed_time_diag seconds.")

# list of phi values used in plotting
unique_phis = p ./ q_list



# =============================== PLOTTING ===============================
start_time_plot = time();



# save plots


# end_time_plot = time();
# elapsed_time_plot = round(end_time_plot - start_time_plot; digits = 3);
# println("All plotting done in $elapsed_time_plot seconds.")
# elapsed_time_all = round(end_time_plot - start_time_init; digits = 3);
# println("Code finished running in $elapsed_time_all seconds. Output files can be found in $plot_save_folder_path.")
