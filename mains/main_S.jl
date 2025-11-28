# calculates the spectrum of a system of a set number of Landau levels (NLL+1)
# in a square cosine periodic potential of given strength U0 [eV] and lattice constant a [m]
# for a given range of rational magnetic flux per unit cell phi = p/q

# arguments can include options for outputting plots and/or data files, performing Wannier analysis,
# simplifying the p/q when possible, varying the number of LLs used at higher fluxes

# using 32-bit floats for memory efficiency and faster diagonalization

#args = ARGS;
args = ["0.021", "5", "8", "120", "0.25", "1", "-p", "-w"]; # for visual studio code testing # "--XBF", "2"

# -OUTPUT FOLDER!-
outfolder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc"
#outfolder = "/users/ivoga/lh/out" # cluster path

outfolder_plots = joinpath(outfolder,"plots/")
outfolder_data = joinpath(outfolder,"data/")
mkpath(outfolder_plots)
mkpath(outfolder_data)


# PACKAGES AND MODULES ===================================================
# packages ------------
using LinearAlgebra
using NPZ
using Plots
using Base.Threads
using Primes
using ProgressMeter
using Printf


# modules -------------
include(joinpath(dirname(@__DIR__),"funcs/aux.jl"))
using .Aux                       # interpret parameters from ARGS; phys consts; set up lists; other auxiliary functions
include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis
include(joinpath(dirname(@__DIR__),"funcs/wannier.jl"))
# modules to be added only if doing wannier analysis or plotting
using .Wannier                      # functions to build and analyse a Wannier plot     
include(joinpath(dirname(@__DIR__),"funcs/plotting.jl"))
using .Plt                          # functions for different kinds of plots   


# PARSE ARGUMENTS AND SET UP LISTS =========================================
parsed = Aux.my_parse(args)
Aux.print_startup_message(parsed)
# bools
plotQ, dataQ, wannierQ, varyNLLQ, plotWQ = parsed.plotQ, parsed.dataQ, parsed.wannierQ, parsed.varyNLLQ, parsed.plotWQ

# ints restricting output to certain LLs or butterflies
XLL, XBF = parsed.XLL, parsed.XBF # 0 by default -- plots full spectrum

# parameters
U0, a_nm, LLmax, q, phi_s, phi_f = parsed.U0, parsed.a, parsed.LLmax, parsed.q, parsed.phi_s, parsed.phi_f
a = a_nm*Float32(1e-9) # convert lattice const to meters

# set up list of fluxes phi = p/q to calculate over
p_list = Aux.get_p_list(q, phi_s, phi_f)
num_fluxes = length(p_list)

# empty arrays to store results and keep track of parameters at each flux after simplifications
phi_list = Vector{Float32}(undef, num_fluxes);  # fluxes
E_list = Vector{Vector{Float32}}(undef, num_fluxes);  # energies
qn_list = Vector{Int64}(undef, num_fluxes);  # simplified qn at each flux
NXY_list = Vector{Int64}(undef, num_fluxes);  # number of XY points at each flux
varyNLLQ ? (NLL_list = Vector{Int64}(undef, num_fluxes)) : (NLL_list = LLmax.*ones(num_fluxes))# number of LLs used at each flux (if varyNLLQ==true)

# number of threads
nthreads = Threads.nthreads()

# MAIN CALCULATION LOOP OVER FLUXES ========================================
println("\nCalculating spectrum over flux range using $nthreads thread(s) for parallel computation.\n")
@showprogress Threads.@threads for j in 1:num_fluxes
    p = p_list[j]
    gcd_pq = gcd(p,q)
    qn = Int64(div(q,gcd_pq)) # simplified qn
    pn = Int64(div(p,gcd_pq))   # simplified pn
    phi = Float32(pn)/Float32(qn) # simplified flux
    phi_list[j] = phi
    qn_list[j] = qn

    energies_phi = Float32[] # energies at this flux, to be pushed to E_list

    NXY = Aux.get_NXY(qn, q) # number of XY-points in each direction
    NXY_list[j] = NXY
    Y_list = collect(range(0f0, stop=Float32(2π), length=(NXY+1)))[1:end-1]
    X_list = Y_list.*Float32(qn)

    #get maximum Landau level index to use
    NLL = varyNLLQ ? Aux.get_NLL_at_flux(LLmax, phi_s, phi) : LLmax
    varyNLLQ && (NLL_list[j] = NLL) 

    for X in X_list
        for Y in Y_list
            H = Hamil.get_full_ham(phi, X, Y, U0, a, pn, NLL)
            evals = eigvals(H)
            append!(energies_phi, Float32.(real.(evals)))
        end
    end
    
    sort!(energies_phi)
    E_list[j] = energies_phi # thread-safe storing of energies

end
# sort all lists by increasing flux
sorted_indices = sortperm(phi_list)
phi_list = phi_list[sorted_indices]
E_list = E_list[sorted_indices]
qn_list = qn_list[sorted_indices]
NXY_list = NXY_list[sorted_indices]
NLL_list = NLL_list[sorted_indices]

# CUT SPECTRUM? ============================================================
out_energies = Aux.cut_spectrum(XLL, XBF, num_fluxes, phi_list, qn_list, NLL_list, NXY_list, E_list)

# OUTPUT PLOTS AND/OR DATA FILES ============================================

# name to attach to filenames, giving the input arguments
lup = varyNLLQ ? "vLL" : "fLL"
param_str = "U0_$(@sprintf("%.5f", U0))_a_$(a_nm)_LL_$(LLmax)_q_$(q)_ps_$(phi_s)_pf_$(phi_f)_XLL_$(XLL)_XBF_$(XBF)_$lup"

if dataQ
    # save energies
    npzwrite(joinpath(outfolder_data, string("energies_",param_str,".npz")), out_energies)
    # save parameters
    npzwrite(joinpath(outfolder_data, string("params_",param_str,".npz")), Dict("phi_list" => phi_list,
                                                           "qn_list" => qn_list,
                                                           "NXY_list" => NXY_list,
                                                           "NLL_list" => NLL_list))
end
if plotQ
    plot_spectrum = Plt.plot_bare_spectrum(out_energies, phi_list)
    Plt.plot_add_LL_guide!(plot_spectrum, phi_s, phi_f, a, LLmax) # saved later in case it needs coloring
end

# WANNIER ANALYSIS (optional) ==============================================
if wannierQ
    wannier_dict = Wannier.mainW(out_energies, phi_list, qn_list, NXY_list)
    # save wannier data
    if dataQ
        npzwrite(joinpath(outfolder_data, "wannier_lines_dict.npz"), wannier_dict)
    end
    if plotWQ
        plot_wannier = Plt.plot_wannier_all(wannier_dict, phi_f, LLmax)
        savefig(plot_wannier, joinpath(outfolder_plots, string("wannier_",param_str,".png")))
    end
    if plotQ
        Plt.color_gaps_eq!(plot_spectrum, wannier_dict, phi_list, 14)        
    end
end



# SAVE FINAL PLOT ==========================================================
if plotQ
    # save spectrum plot
    savefig(plot_spectrum, joinpath(outfolder_plots, string("spectrum_",param_str,".png")))
end



