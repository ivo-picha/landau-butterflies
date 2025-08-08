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

U0 = 0.15 # potential strenght in eV
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
    append!(energies,evals)
    ktot = norm(kpoint)*sqrt(2)*a/(2π)
    append!(ktots,[ktot for _ in eachindex(evals)])
end

pbs = scatter(ktots,energies, msw=0, color = :red, label = "", ms = 0.7, ylims = (-0.1,0.65), framestyle = :box, xlabel = "|k|", ylabel = "ε")
title!(pbs, "U₀ = $U0 eV, a = $(a*1e9) nm")
#save plot
plot_name = "BS_ZF_U$U0.png"
savefig(pbs, joinpath(plot_folder,plot_name))
