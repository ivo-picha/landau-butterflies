# calculates the spectrum of a system of a set number of Landau levels (NLL+1)
# in a square cosine periodic potential of given strength U0 [eV] and lattice constant a [m]
# for a given rational magnetic flux per unit cell phi = p/q

# calculates the density of states OF THE LOWEST BAND and outputs a plot and or a data file
using Plots
using NPZ
using LinearAlgebra
using ProgressMeter
using Base.Threads
using KernelDensity
using Measures

include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis

# output folder
outfolder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc"
outfolder_plots = joinpath(outfolder,"plots_DOS/")
outfolder_data = joinpath(outfolder,"data_DOS/")
mkpath(outfolder_plots)
mkpath(outfolder_data)

# PARAMETERS ============================================================

p = 1
q = 2
U0 = 0.05f0          # potential strength [eV]
a_nm = 5.0f0         # lattice constant [nm]
NLL = Int64(q/p * 30)              # number of Landau levels to include + 1

a = a_nm*1f-9       # convert lattice const to meters
NXY = 129          # number of XY in each direction in MBZ


if gcd(p,q) != 1
    p = div(p, gcd(p,q))
    q = div(q, gcd(p,q))
end

phi = Float32(p/q)

X_list = range(0f0, Float32(2π*q), NXY+1)[1:end-1]
Y_list = range(0f0, Float32(2π), NXY+1)[1:end-1]
zipped_list = reshape(collect(Iterators.product(X_list, Y_list)),:)


energies = Array{Float32}(undef, q, length(zipped_list))
nt = Threads.nthreads()
@showprogress Threads.@threads for i in 1:NXY^2
    X = zipped_list[i][1]
    Y = zipped_list[i][2]
    H = Hamil.get_full_ham(phi, X, Y, U0, a, p, NLL)
    energies[:,i] = eigvals(H)[1:q]
end

energies = sort(reshape(energies, :,))
eta = (maximum(energies)-minimum(energies))/100         # broadening for DOS [eV]
kd = kde(energies, bandwidth = eta)

# PLOTTING AND OUTPUT ====================================================
plt = plot(
    kd.x,
    kd.density,
    linewidth = 2,
    xlabel = "Energy",
    ylabel = "Density of States",
    yticks = false,
    legend = false,
    framestyle = :box,
    title = "DOS for U0=$(U0)eV, a=$(a_nm)nm, φ=$(p)/$(q), NLL=$(NLL)",
    size = (600,400),
    margin = 5mm,
)
png_filename = joinpath(outfolder_plots, "DOS_U0=$(U0)_a=$(a_nm)_phi=$(p)_$(q)_NLL=$(NLL).png")
savefig(plt, png_filename)

# save data
data_filename = joinpath(outfolder_data, "DOS_U0=$(U0)_a=$(a_nm)_phi=$(p)_$(q)_NLL=$(NLL).npy")
npzwrite(data_filename, energies)
