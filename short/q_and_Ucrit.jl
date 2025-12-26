# calculates the spectrum of a system of a set number of Landau levels (NLL+1)
# in a square cosine periodic potential of given strength U0 [eV] and lattice constant a [m]
# for a given range of rational magnetic flux per unit cell phi = 1/qmax

# finds critical U above which a Hofstadter butterfly gap opens

# using 32-bit floats for memory efficiency and faster diagonalization

args = ARGS;
#args = ["0.01", "5", "4"]; # for visual studio code testing 

# -OUTPUT FOLDER!-
#outfolder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc"
outfolder = "/users/ivoga/lh/out" # cluster path

outfolder_plots = joinpath(outfolder,"plots/")
outfolder_data = joinpath(outfolder,"data/")
mkpath(outfolder_plots)
mkpath(outfolder_data)


# PACKAGES AND MODULES ===================================================
# packages ------------
using LinearAlgebra
using Plots
using Base.Threads
using ProgressMeter
using Printf
using Measures
using NPZ

include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis

# parameters
U0, a_nm, q_max = parse.(Float32, args)
q_max = Int64(q_max)
a = a_nm*Float32(1e-9) # convert lattice const to meters

# set up list of fluxes phi = p/q to calculate over
q_list = range(1, q_max)
num_fluxes = length(q_list)

# empty arrays to store results and keep track of parameters at each flux after simplifications
qn_list = Vector{Int64}(undef, num_fluxes);  # denominators
E_list = Vector{Vector{Float32}}(undef, num_fluxes);  # energies
NXY = 32  # number of X points and Y points at each flux

# number of threads
nthreads = Threads.nthreads()

# MAIN CALCULATION LOOP OVER FLUXES ========================================
println("\nCalculating spectrum over flux range using $nthreads thread(s) for parallel computation.\n")
@showprogress Threads.@threads for j in 1:num_fluxes
    qn = q_list[j]
    qn_list[j] = qn
    phi = Float32(1f0)/Float32(qn)

    energies_phi = Float32[] # energies at this flux, to be pushed to E_list

    Y_list = collect(range(0f0, stop=Float32(2π), length=(NXY+1)))[1:end-1]
    X_list = Y_list.*Float32(qn)

    #get maximum Landau level index to use
    NLL = Int64(clamp(10*qn, 0,500))

    for X in X_list
        for Y in Y_list
            H = Hamil.get_full_ham(phi, X, Y, U0, a, 1, NLL)
            evals = eigvals(H)
            append!(energies_phi, Float32.(real.(evals)))
        end
    end
    
    sort!(energies_phi)
    E_list[j] = energies_phi[1:Int(2*qn*NXY^2)] # thread-safe storing of energies

end
# sort all lists by increasing q
sorted_indices = sortperm(qn_list)
qn_list = qn_list[sorted_indices]
E_list = E_list[sorted_indices]

#PLOTTING ==============================================================
spectrum_bare_options = (
    markersize = 2,
    color = :black,
    label = "",
    xlabel = "q",
    ylabel = "E [meV]",
    framestyle = :box,
    size = (1200,800),
    tickfontsize = 16,
    guidefontsize = 18,
    margin = 9mm
)
plot1 = Plots.Plot()
for (j,ens) in enumerate(E_list)
    xs = [qn_list[j] for i = 1:length(ens)]
    Plots.scatter!(plot1, xs, ens; spectrum_bare_options...)
end
Plots.title!(plot1, "U0 = $(U0) eV, a = $(a_nm) nm")
miny, maxy = minimum(vcat(E_list...)), maximum(vcat(E_list...))
Plots.ylims!(plot1, miny - 0.05f0*(maxy - miny), maxy + 0.05f0*(maxy - miny))

png_path = joinpath(outfolder_plots, "q_spectrum_U0_$(U0)_a_$(a_nm)_1oqlim_$(q_max).png")
println("Saving spectrum plot to $png_path")
Plots.savefig(plot1, png_path)

# SAVING DATA ==============================================================
data_path = joinpath(outfolder_data, "q_spectrum_U0_$(U0)_a_$(a_nm)_1oqlim_$(q_max).npz")
println("Saving spectrum data to $data_path")
# NPZ cannot write a Vector{Vector{Float32}} directly because it's a jagged array
# Save `qn_list` and each inner energy vector as separate arrays inside the .npz
data_dict = Dict{String,Any}()
data_dict["qn_list"] = qn_list
for (j, ens) in enumerate(E_list)
    data_dict[string("E_", j)] = ens
end
npzwrite(data_path, data_dict)
