# code, which calculates the spectrum for a given a, U₀, and flux (and # of LLs)
# and plots the bands in k_y and Y quantum number space, akin to a bandstructure

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
using FileIO

plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/bands"
#plot_save_folder_path = "/users/ivoga/lh/plts/spectra"

#args = ARGS
args = ["[2, 2, 0.03, 50, 14]"]

# get parameters from ARGS
p, q, U0, a_in_angstr, NLL = Params.parse_arguments_bands(args)
if gcd(p,q) != 1
    println("p/q can be simplified. not doing it for u tho")
end
a = a_in_angstr * 1e-10                     # lattice const in meters
phi = p/q
xi0 = sqrt(2π / phi)

# get lists of ky0 and Y values to iterate over
Nky = 64;                                    # number of ky* points; independent calculations; variation on scale of U0
ky_list = Params.get_ky_list(a, Nky+1)
NY = Nky;
Y_list = Params.get_Y_list(NY+1)

# empty lists to store coordinates for plot
nbands = p*(NLL+1)
band_energies = [zeros(Float64, Nky, NY) for _ in 1:nbands]

# message for time it took to initialise
end_time_init = time();
elapsed_time_init = round(end_time_init - start_time_init; digits = 2);
println("Initialisation done in $elapsed_time_init seconds.")

println("Calculating band structure...")

# iterate over lists and diagonalise hamiltonians
start_time_diag = time();                   # set up a clock to monitor elapsed time

@showprogress for (i,ky) in enumerate(ky_list)
    for (j,Y) in enumerate(Y_list)

        H = Hamil.get_full_ham(xi0, ky, Y, U0, a, p, NLL)
        evalsH = eigvals(H)

        for (n,ev) in enumerate(evalsH)
            band_energies[n][i,j] = ev              
        end

    end
end

end_time_diag = time();
elapsed_time_diag = round(end_time_diag - start_time_diag; digits = 2);
println("Band structure has been calculated in $elapsed_time_diag seconds. \nPlotting...")



# =============================== PLOTTING ===============================

xgrid = collect(ky_list).*a # ky*a
ygrid = collect(Y_list) # Y
plotBS = Plt.plot_band_structure(xgrid, ygrid, band_energies, 1)
#last argument tells how many of the bottom bands to plot; leave empty or 0 for all bands

FileIO.save(joinpath(plot_save_folder_path, "bs_p$p-q$q-U$U0-a$a_in_angstr-NLL$NLL-Nk$Nky.png"), plotBS)

end_time_plot = time()
elapsed_time_all = round(end_time_plot - start_time_init; digits = 2);
println("Code finished running in $elapsed_time_all seconds. Output files can be found in $plot_save_folder_path.")
