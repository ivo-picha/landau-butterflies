# code which takes a U₀ and an a, and for a set of fluxes computes the DOS per energy up to
# some specified filling (set to 1 automatically). outputs a 3D lineplot of DOS lines for each flux

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
using FileIO

#plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/SCWC/"
plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/DOS"

#args = ARGS
args = ["[1.5, 2.0, 0.05, 50, 51, 4, 1.0]"]

# get parameters from ARGS
startphi, endphi, U0, a_in_angstr, q, Nmin, np = Params.parse_arguments_DOS(args)
a = a_in_angstr * 1e-10                     # lattice const in meters

# get lists of q, ky0 and Y values to iterate over
Npp = 15 # number of desired lines in phi in the 3d plot
p_list = unique(Int.(round.(collect(range(q*startphi,q*endphi, Npp)))))
Nky = 65;                                    # number of ky* points; independent calculations; variation on scale of U0
ky_list = Params.get_ky_list(a, Nky)
NY = Nky;
Y_list = Params.get_Y_list(NY)


# empty lists to store coordinates for plot
DOSs = Tuple{Vector{Float64},Vector{Float64}}[];                       # DOS at each flux (bin_centers, norm_frequencies)

# set up an energy limit
Emax = Hamil.E_LL(Nmin,sqrt(2π/endphi), a) + U0

# message for time it took to initialise
end_time_init = time();
elapsed_time_init = round(end_time_init - start_time_init; digits = 3);
println("Initialisation done in $elapsed_time_init seconds.")
println("Calculating DOS up to filling $np from flux $startphi to $endphi.")
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

    # get DOS plot vectors at flux
    bins_phi, freq_phi = Wannier.get_DOS_upto_np(energies_at_phi, phi, NLL, np)
    push!(DOSs, (bins_phi, freq_phi))
    

end
end_time_diag = time();
elapsed_time_diag = round(end_time_diag - start_time_diag; digits = 3);
println("Spectrum and DOS have been calculated in $elapsed_time_diag seconds.")

# list of phi values used in plotting
unique_phis = p_list ./ q


# save plots
# spectrum_plot_name = string("SCWC_fixE_vp_S_U$U0-a$a_in_angstr-q$q-phi$startphi-phf$endphi")
# spectrum_plot_path = string(joinpath(plot_save_folder_path, spectrum_plot_name), ".png")
# savefig(plot_s, spectrum_plot_path)


pltDOS = Plt.plot_3D_DOS(DOSs, unique_phis)

FileIO.save(joinpath(plot_save_folder_path, "dos_U$U0-a$a_in_angstr-n$np-Nk$Nky.svg"), pltDOS)

end_time = time();
elapsed_time_all = round(end_time - start_time_init; digits = 3);
println("Code finished running in $elapsed_time_all seconds. Output files can be found in $plot_save_folder_path.")

