# find the spectrum of Landau levels in a cosine potential in 2D as a function of flux
# set a max E and calculate different number of landau levels for each flux; vary p to simplest coprime with q
# only plot 

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
using NPZ

#plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/S_LB/"
plot_save_folder_path = "/users/ivoga/lh/plts/spectra"
#data_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/local/LB_spectra/"
data_save_folder_path = "/users/ivoga/lh/data_LB"

args = ARGS
#args = ["[0.25, 0.5, 0.01, 50, 48, 1, 0.05]"]

# get parameters from ARGS
startphi, endphi, U0, a_in_angstr, q, Nmin, gap_factor = Params.parse_arguments(args)
a = a_in_angstr * 1e-10                     # lattice const in meters

# get lists of q, ky0 and Y values to iterate over
p_list = unique(Int.(collect(range(round(q*startphi),round(q*endphi)))))
Nky = 16;                                    # number of ky* points; independent calculations; variation on scale of U0
ky_list = Params.get_ky_list(a, Nky)
NY = Nky;
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

# set up an energy limit
Emax = Hamil.E_LL(Nmin,sqrt(2π/endphi), a) + U0

# message for time it took to initialise
end_time_init = time();
elapsed_time_init = round(end_time_init - start_time_init; digits = 3);
println("Initialisation done in $elapsed_time_init seconds.")

Params.startmessage_S_fixE(startphi, endphi, U0, a_in_angstr, q, Nmin)
println(" p/q is simplified to an irrational fraction.")

# iterate over lists and diagonalise hamiltonians
start_time_diag = time();                   # set up a clock to monitor elapsed time
nt = nthreads() # number of threads
@showprogress for p in p_list

    gcdpq = gcd(p,q)
    pn = Int(p/gcdpq)
    qn = Int(q/gcdpq)

    phi = pn/qn                               # unit flux per unit cell
    xi0 = sqrt(2π / phi)
    energies_at_phi = Float64[];            # to be appended to global energies list

    # find the max NLL to be plotted at that energy;
    # never go past 30
    local NLL
    NLL::Int = 0;
    while Hamil.E_LL(NLL,xi0,a) < Emax && NLL < 31
        NLL += 1
    end

    # set up list for each thread
    tlists = [Vector{Float64}() for _ in 1:nt]

    @threads for ky in ky_list
        tid = threadid()

        for Y in Y_list

        H = Hamil.get_full_ham(xi0, ky, Y, U0, a, pn, NLL)
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
unique_phis = p_list ./ q

# vector which contains full information for each point of Wannier plot as tuple (phi, pdensity, gap_size, gap_lower_energy)
wannier_points = collect(zip(phis_w, Float64.(particle_densities_global), Float64.(gaps_global), Float64.(energies_lower_gaps_global)))
# dictionary for unique lines in Wannier plot
lines_dict = Dict{Tuple{Float64, Float64}, Vector{NTuple{4, Float64}}}()

# update dictionary with all available lines
Wannier.identify_lines(lines_dict, unique_phis, wannier_points, phis_w)

# merge lines that are similar enough
merged_dict = Wannier.merge_round_keys(lines_dict)

# =============================== PLOTTING ===============================
start_time_plot = time();
plots_title = string("U₀=$U0 eV,  a=$a_in_angstr Å")

# maximum energy of the lowest band (defined by 1 particle per unit cell); take only points in LB for simplicity
Emax_LB = maximum(t[4] for t in merged_dict[(0.0,1.0)])
enmask = energies .< Emax_LB
phis_LB = phis[enmask]
energies_LB = energies[enmask]

# save data so it can be accessed later; npz format, readable by python as well
npzwrite(joinpath(data_save_folder_path, "LB_S_U$U0-a$a_in_angstr-q$q-phi$startphi-phf$endphi.npz"),
        Dict("x" => phis_LB, "y" => energies_LB))

# plot only the spectrum
plot_s = Plt.plot_spectrum_bare(phis, energies, plots_title, (0.0, endphi), (minimum(energies_LB), Emax_LB))



# save plots
spectrum_plot_name = string("SCWC_fixE_vp_LB_S_U$U0-a$a_in_angstr-q$q-phi$startphi-phf$endphi")
spectrum_plot_path = string(joinpath(plot_save_folder_path, spectrum_plot_name), ".png")
savefig(plot_s, spectrum_plot_path)


end_time_plot = time();
elapsed_time_plot = round(end_time_plot - start_time_plot; digits = 3);
println("All plotting done in $elapsed_time_plot seconds.")
elapsed_time_all = round(end_time_plot - start_time_init; digits = 3);
println("Code finished running in $elapsed_time_all seconds. Output files can be found in $plot_save_folder_path.")

