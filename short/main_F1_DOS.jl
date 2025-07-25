# find the spectrum and DOS of the lowest band for a given flux

include(joinpath(dirname(@__DIR__),"mods/params.jl"))
using .Params                       # interpret parameters from ARGS; phys consts; set up lists
include(joinpath(dirname(@__DIR__),"mods/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis                   # wavefunctions and operations for them + density plotting
include(joinpath(dirname(@__DIR__),"mods/plotting.jl"))
using .Plt                          # functions for different kinds of plots      
include(joinpath(dirname(@__DIR__),"mods/densities.jl"))
using .Dens


using ProgressMeter
using LinearAlgebra
using Plots


plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/random/"


U0 = 0.04
p = 3
q = 1
a_in_angstr = 50
a = a_in_angstr * 1e-10                     # lattice const in meters

NLL = 40

phi = p/q                                   # unit flux per unit cell
xi0 = sqrt(2π / phi)

# get lists of ky0 and Y values to iterate over
Nky = 32;                                   # number of ky* points; independent calculations; variation on scale of U0
ky_list = Params.get_ky_list(a, Nky)
NY = 32;
Y_list = Params.get_Y_list(NY)


# iterate over lists and diagonalise hamiltonians
start_time_diag = time();                   # set up a clock to monitor elapsed time

energies = Float64[];

@showprogress for Y in Y_list
    for ky0 in ky_list

        H = Hamil.get_full_ham(xi0, ky0, Y, U0, a, p, NLL)
        evalsH = eigvals(H)
        append!(energies, evalsH)
    end
end
sort!(energies)
energies = energies[1:(q*(Nky-1)*(NY-1))]

# gaps = diff(energies)
# gappos = findfirst(x-> x>U0/2, gaps)
# energies_LB = energies[1:gappos]
# sort!(energies)

# gaps = diff(energies)
# gappos = findfirst(x-> x>U0/2, gaps)
# energies_LB = energies[1:gappos]

# make a histogram output out of intput of energy list, for DOS
function histogram_data(data::Vector{Float64}, N::Int, norm::Float64=1.0)
    # Compute min and max of the data
    min_val = data[1]
    max_val = data[end]
    
    # Compute bin edges and centers
    edges = range(min_val, max_val; length=N+1)
    bin_centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])

    # Initialize frequency vector
    counts = zeros(Int, N)

    # Bin the data manually
    for x in data
        if x == max_val
            # put max value in last bin
            counts[end] += 1
        elseif x == min_val
            counts[1] += 1
        else
            bin_index = searchsortedfirst(edges, x)
            counts[bin_index - 1] += 1
        end
    end

    # Normalize frequencies to sum to np
    total = sum(counts)
    frequencies = counts .* (norm / total)

    return (bin_centers, frequencies)
end

bins, freqs = histogram_data(energies, 25)

p2 = Plots.plot(bins,freqs,
    label = "", framestyle=:box, xlabel="E [eV]", ylabel="DOS", color = :blue,
    ylims=(0,1.1*maximum(freqs)), title = "ϕ = $p/$q, U₀ = $U0 eV, a = $(a*1e10) Å")

#save plot
plot_name = "DOS_F$phi-U$U0-N$NLL.png"
savefig(p2, joinpath(plot_save_folder_path,plot_name))
