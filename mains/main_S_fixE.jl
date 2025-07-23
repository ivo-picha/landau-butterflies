# find the spectrum of Landau levels in a cosine potential in 2D as a function of flux
# set a max E and calculate different number of landau levels for each flux 

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

plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/SCWC/"
#plot_save_folder_path = "/users/ivoga/lh/plts/spectra"

#args = ARGS
args = ["[0.1, 2.0, 0.01 , 50, 51, 3, 0.0]"]

# get parameters from ARGS
startphi, endphi, U0, a_in_angstr, p, Nmin = Params.parse_arguments(args)[1:end-1]
a = a_in_angstr * 1e-10                     # lattice const in meters

# get lists of q, ky0 and Y values to iterate over
q_list = Params.get_q_list_red(startphi, endphi, p)
Nky = 10;                                    # number of ky* points; independent calculations; variation on scale of U0
ky_list = Params.get_ky_list(a, Nky)
NY = 10;
Y_list = Params.get_Y_list(NY)


# empty lists to store coordinates for plot
energies = Float64[];                       # all energies 
phis = Float64[];                           # repeating phis, for easier plotting later on


# set up an energy limit
Emax = Hamil.E_LL(Nmin,sqrt(2π/endphi), a) + U0

# message for time it took to initialise
end_time_init = time();
elapsed_time_init = round(end_time_init - start_time_init; digits = 3);
println("Initialisation done in $elapsed_time_init seconds.")

Params.startmessage_S_fixE(startphi, endphi, U0, a_in_angstr, p, Nmin)

# iterate over lists and diagonalise hamiltonians
start_time_diag = time();                   # set up a clock to monitor elapsed time
@showprogress for q in q_list

    phi = p/q                               # unit flux per unit cell
    xi0 = sqrt(2π / phi)
    energies_at_phi = Float64[];            # to be appended to global energies list

    # find the max NLL to be plotted at that energy;
    # never go past 30
    NLL::Int = 0;
    while Hamil.E_LL(NLL,xi0,a) < Emax && NLL < 31
        NLL += 1
    end


    for ky in ky_list
        for Y in Y_list

        H = Hamil.get_full_ham(xi0, ky, Y, U0, a, p, NLL)
        evalsH = eigvals(H)
        energies_at_phi = [energies_at_phi; evalsH] #add eigenvalues to list of energies

        end
    end

    global energies
    energies = [energies; energies_at_phi]
    global phis
    phis = [phis; [phi for j = 1:length(energies_at_phi)]] #add phi values to list of phis

end
end_time_diag = time();
elapsed_time_diag = round(end_time_diag - start_time_diag; digits = 3);
println("Spectrum has been calculated in $elapsed_time_diag seconds.")




# =============================== PLOTTING ===============================
start_time_plot = time();
plots_title = string("U₀=$U0 eV,  a=$a_in_angstr Å")


# plot only the spectrum
plot_s = Plt.plot_spectrum_bare(phis, energies, plots_title, (0.0, endphi), (minimum(energies), Emax))

# add guiding lines
Plt.plot_add_LL_guide!(plot_s, startphi, endphi, a, 30)


# save plots
spectrum_plot_name = string("S_fixE_U$U0-a$a_in_angstr-p$p-phi$startphi-phf$endphi")
spectrum_plot_path = string(joinpath(plot_save_folder_path, spectrum_plot_name), ".png")
savefig(plot_s, spectrum_plot_path)


end_time_plot = time();
elapsed_time_plot = round(end_time_plot - start_time_plot; digits = 3);
println("All plotting done in $elapsed_time_plot seconds.")
elapsed_time_all = round(end_time_plot - start_time_init; digits = 3);
println("Code finished running in $elapsed_time_all seconds. Output files can be found in $plot_save_folder_path.")

