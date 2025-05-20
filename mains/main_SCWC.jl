# find the spectrum of Landau levels in a cosine potential in 2D as a function of flux
# find the corresponding Wannier plot and color gaps in the spectrum according to Chern number

start_time_init = time();

include(joinpath(dirname(@__DIR__),"mods/params.jl"))
using .Params                       # interpret parameters from ARGS; phys consts; set up lists
include(joinpath(dirname(@__DIR__),"mods/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis
include(joinpath(dirname(@__DIR__),"mods/wannier.jl"))
using .Wannier                      # functions to build and analyse a Wannier plot     
include(joinpath(dirname(@__DIR__),"mods/plotting.jl"))
using .Plt                          # functions for different kinds of plots      

using ProgressMeter
using LinearAlgebra
using Plots
using Base.Threads

plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/SCWC/"
#plot_save_folder_path = "/users/ivoga/lh/plts/spectra"

#args = ARGS
args = ["[0.4, 0.5, 0.01 , 50, 21, 30, 0.2]"]

# get parameters from ARGS
startphi, endphi, U0, a_in_angstr, p, NLL, gap_factor = Params.parse_arguments(args)
a = a_in_angstr * 1e-10                     # lattice const in meters

# get lists of q, ky0 and Y values to iterate over
q_list = Params.get_q_list(startphi, endphi, p)
Nky = 10;                                    # number of ky* points; independent calculations; variation on scale of U0
ky_list = Params.get_ky_list(a, Nky)
NY = 10;
Y_list = Params.get_Y_list(NY)

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

Params.startmessage_SCWC(startphi, endphi, U0, a_in_angstr, p, NLL)
# message for size of calculation
Params.print_size_message(q_list, p, Nky, NY, NLL)

# iterate over lists and diagonalise hamiltonians
start_time_diag = time();                   # set up a clock to monitor elapsed time
nt = nthreads() # number of threads
@showprogress for q in q_list

    phi = p/q                               # unit flux per unit cell
    xi0 = sqrt(2π / phi)
    energies_at_phi = Float64[];            # to be appended to global energies list

    # set up list for each thread
    tlists = [Vector{Float64}() for _ in 1:nt]

    @threads for ky in ky_list
        tid = threadid()

        for Y in Y_list

        H = Hamil.get_full_ham(xi0, ky, Y, U0, a, p, NLL)
        evalsH = eigvals(H)
        append!(tlists[tid], evalsH) #add eigenvalues to list of energies

        end
    end

    # combine buffer lists
    energies_at_phi = reduce(vcat, tlists)

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

# vector which contains full information for each point of Wannier plot as tuple (phi, pdensity, gap_size, gap_lower_energy)
wannier_points = collect(zip(phis_w, Float64.(particle_densities_global), Float64.(gaps_global), Float64.(energies_lower_gaps_global)))
# dictionary for unique lines in Wannier plot
lines_dict = Dict{Tuple{Float64, Float64}, Vector{NTuple{4, Float64}}}()

# update dictionary with all available lines
Wannier.identify_lines(lines_dict, unique_phis, wannier_points, phis_w)

merged_dict = Wannier.merge_round_keys(lines_dict)

# =============================== PLOTTING ===============================
start_time_plot = time();
plots_title = string("U₀=$U0 eV,  a=$a_in_angstr Å,  Nₗₗ=$NLL")

# plot colored Wannier plot
plot_w = Plt.plot_wannier_all(wannier_points, gaps_global, merged_dict, endphi, NLL, plots_title)
title!(plot_w, plots_title)

# plot only the spectrum
plot_s = Plt.plot_spectrum_bare(phis, energies, plots_title)
title!(plot_s, plots_title)

# plot colors in the gaps of the spectrum
Plt.color_gaps!(plot_s, merged_dict, unique_phis, NLL)

# add guiding lines
Plt.plot_add_LL_guide!(plot_s, startphi, endphi, a, NLL)


# save plots
spectrum_plot_name = string("SCWC_S_N$NLL-U$U0-a$a_in_angstr-p$p-phi$startphi-phf$endphi")
spectrum_plot_path = string(joinpath(plot_save_folder_path, spectrum_plot_name), ".png")
savefig(plot_s, spectrum_plot_path)
wannier_plot_name = string("SCWC_W_N$NLL-U$U0-a$a_in_angstr-p$p-phi$startphi-phf$endphi")
wannier_plot_path = string(joinpath(plot_save_folder_path, wannier_plot_name), ".png")
savefig(plot_w, wannier_plot_path)


end_time_plot = time();
elapsed_time_plot = round(end_time_plot - start_time_plot; digits = 3);
println("All plotting done in $elapsed_time_plot seconds.")
elapsed_time_all = round(end_time_plot - start_time_init; digits = 3);
println("Code finished running in $elapsed_time_all seconds. Output files can be found in $plot_save_folder_path.")


# for (key,value) in merged_dict
#     if 0.98 < key[1] < 1.02 && abs(key[2])<0.001
#         println("$key : $value")
#     end
        
# end

# length(merged_dict[(1.0,0.0)])
# length(lines_dict[(1.0,0.0)])
# length(unique_phis)