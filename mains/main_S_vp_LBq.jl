# find the spectrum of Landau levels in a cosine potential in 2D as a function of flux
# set a max E and calculate different number of landau levels for each flux; vary p to simplest coprime with q
# only plot the lowest q bands (lowest band -- hofstadter)

start_time_init = time();

include(joinpath(dirname(@__DIR__),"mods/params.jl"))
using .Params                       # interpret parameters from ARGS; phys consts; set up lists
include(joinpath(dirname(@__DIR__),"mods/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis
include(joinpath(dirname(@__DIR__),"mods/plotting.jl"))
using .Plt                          # functions for different kinds of plots      

using ProgressMeter
using LinearAlgebra
using Plots
# using Base.Threads
using NPZ

#plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/S_LB/"
plot_save_folder_path = "/users/ivoga/lh/plts/spectra"
#data_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/local/LB_spectra/"
data_save_folder_path = "/users/ivoga/lh/data_LB"

args = ARGS
#args = ["[0.25, 0.5, 0.05, 50, 24, 40]"]

# get parameters from ARGS
startphi, endphi, U0, a_in_angstr, q, NLL = Params.parse_arguments_vp_LBq(args)
a = a_in_angstr * 1e-10                     # lattice const in meters

# get lists of q, ky0 and Y values to iterate over
p_list = unique(Int.(collect(range(round(q*startphi),round(q*endphi)))))
Nky = 16;                                    # number of ky* points; independent calculations; variation on scale of U0
ky_list = Params.get_ky_list(a, Nky)
NY = Nky;
Y_list = Params.get_Y_list(NY)


# empty lists to store coordinates for plot
energies_LB = Float64[];                       # all energies 
phis = Float64[];                           # repeating phis, for easier plotting later on


# message for time it took to initialise
end_time_init = time();
elapsed_time_init = round(end_time_init - start_time_init; digits = 3);
println("Initialisation done in $elapsed_time_init seconds.")

Params.startmessage_S_vp(startphi, endphi, U0, a_in_angstr, q, NLL)
println(" p/q is simplified to an irrational fraction.")

# iterate over lists and diagonalise hamiltonians
start_time_diag = time();                   # set up a clock to monitor elapsed time
# nt = nthreads() # number of threads
@showprogress for p in p_list

    gcdpq = gcd(p,q)
    pn = Int(p/gcdpq)
    qn = Int(q/gcdpq)

    phi = pn/qn                               # unit flux per unit cell
    xi0 = sqrt(2π / phi)
    energies_at_phi = Float64[];            # to be appended to global energies list

    # # set up list for each thread
    # tlists = [Vector{Float64}() for _ in 1:nt]

    # @threads for ky in ky_list
    #     tid = threadid()

    #     for Y in Y_list

    #     H = Hamil.get_full_ham(xi0, ky, Y, U0, a, pn, NLL)
    #     evalsH = eigvals(H)
    #     append!(tlists[tid], evalsH) #add eigenvalues to list of energies

    #     end
    # end

    # # combine buffer lists
    # energies_at_phi = reduce(vcat, tlists)

    for ky in ky_list
        for Y in Y_list

        H = Hamil.get_full_ham(xi0, ky, Y, U0, a, pn, NLL)
        evalsH = eigvals(H)
        append!(energies_at_phi, evalsH) #add eigenvalues to list of energies

        end
    end

    sort!(energies_at_phi)
    energies_at_phi_LB = energies_at_phi[1:(qn*(Nky-1)*(NY-1))]

    global energies_LB
    append!(energies_LB, energies_at_phi_LB)
    global phis
    phis = [phis; [phi for j = 1:length(energies_at_phi_LB)]] #add phi values to list of phis

end
end_time_diag = time();
elapsed_time_diag = round(end_time_diag - start_time_diag; digits = 3);
println("Spectrum has been calculated in $elapsed_time_diag seconds.")


# =============================== PLOTTING ===============================
start_time_plot = time();
plots_title = string("U₀=$U0 eV,  a=$a_in_angstr Å")

# save data so it can be accessed later; npz format, readable by python as well
npzwrite(joinpath(data_save_folder_path, "LB_S_U$U0-a$a_in_angstr-q$q-phi$startphi-phf$endphi-N$NLL.npz"),
        Dict("x" => phis, "y" => energies_LB))

# # plot only the spectrum
# plot_s = Plt.plot_spectrum_bare(phis, energies_LB, plots_title, (0.0, endphi), (minimum(energies_LB), maximum(energies_LB)));


# # save plots
# spectrum_plot_name = string("S_vp_LBq_U$U0-a$a_in_angstr-q$q-phi$startphi-phf$endphi-N$NLL")
# spectrum_plot_path = string(joinpath(plot_save_folder_path, spectrum_plot_name), ".png")
# savefig(plot_s, spectrum_plot_path)


end_time_plot = time();
elapsed_time_plot = round(end_time_plot - start_time_plot; digits = 3);
println("All plotting done in $elapsed_time_plot seconds.")
elapsed_time_all = round(end_time_plot - start_time_init; digits = 3);
println("Code finished running in $elapsed_time_all seconds. Output files can be found in $plot_save_folder_path.")

