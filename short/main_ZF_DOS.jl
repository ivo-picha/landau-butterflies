# spectrum and DOS for 2d electron gas in cos potential
using LinearAlgebra
using ProgressMeter
using Plots
using Statistics


plot_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/random"
#plot_folder = "/users/ivoga/lh-short/plts"

const ħ = 6.62607015e-34/(2π);  # Planck constant [J s]
const e = 1.602176634e-19;      # elementary charge [C]
const m_e = 9.1093837139e-31;   # electron mass [kg];

U0 = 0.05 # potential strenght in eV
a = 5e-9; # lattice constant
m = m_e; # electron/particle mass

G = 2π/a # recip scat vec
# nu = 1.0; # filling
# EF = (ħ^2 /e)* 2π*(nu/a^2)/(2*m) 

Nk = 100; # sqrt of number of momentum states in BZ
NBZ = 5; # number of BZs in each dimension / 2

BZ_centers = reshape(collect(Base.product(G.*(-NBZ:NBZ), G.*(-NBZ:NBZ))), :)
BZ_kpoints = reshape(collect(Base.product(range(-G/2,G/2,Nk), range(-G/2,G/2,Nk))), :)

ϵ = G/1000 # numerical error tolerance

function get_Hk(BZ_kpoint::Tuple{Float64,Float64}, BZ_centers::Vector{Tuple{Float64,Float64}}, m::Float64, ϵ::Float64)
    # diagonal elements
    Hk = (ħ^2 /e).*diagm([((BZ_kpoint[1] + kc[1])^2 + (BZ_kpoint[2] + kc[2])^2)/(2*m) for kc in BZ_centers])
    # off-diagonal elements
    for i in eachindex(BZ_centers)
        for j in eachindex(BZ_centers)
            if abs(norm(BZ_centers[i] .- BZ_centers[j]) - G) < ϵ
                Hk[i,j] = U0/2
            end
        end
    end
    return Hermitian(Hk)
end

println("Calculating spectrum and eigenstates...")

energies = Float64[];
ktots = Float64[];
@showprogress for kpoint in BZ_kpoints
    Hk = get_Hk(kpoint, BZ_centers, m, ϵ)
    evals = eigvals(Hk)
    for i = eachindex(evals)
        append!(energies, evals)
    end
end

sort!(energies)

gaps = diff(energies)
gappos = findfirst(x-> x>U0/2, gaps)
energies[gappos-1]
energies_LB = energies[1:gappos]

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

bins, freqs = histogram_data(energies_LB, 25)

p1 = Plots.plot(bins,freqs,
    label = "", framestyle=:box, xlabel="E [eV]", ylabel="DOS", color = :red,
    ylims=(0,maximum(freqs)), title = "U₀ = $U0 eV, a = $(a*1e10) Å")

#save plot
# plot_name = "DOS_ZF_U$U0.png"
# savefig(p1, joinpath(plot_folder,plot_name))
