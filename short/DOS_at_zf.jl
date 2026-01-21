# spectrum and DOS for 2d noninteracting electron gas in cos potential; no magnetic field
using LinearAlgebra
using ProgressMeter
using Plots
using Statistics
using KernelDensity
using Measures
using NPZ


# output folder
outfolder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc"
outfolder_plots = joinpath(outfolder,"plots_DOS/")
outfolder_data = joinpath(outfolder,"data_DOS/")
mkpath(outfolder_plots)
mkpath(outfolder_data)

const ħ = 6.62607015f-34/(2π);  # Planck constant [J s]
const e = 1.602176634f-19;      # elementary charge [C]
const m_e = 9.1093837139f-31   # electron mass [kg];

U0 = 0.05f0 # potential strenght in eV
a = 5f-9; # lattice constant
m = m_e; # electron/particle mass

G = Float32(2π/a) # recip scat vec
# nu = 1.0; # filling
# EF = (ħ^2 /e)* 2π*(nu/a^2)/(2*m) 

Nk = 100; # sqrt of number of momentum states in BZ
NBZ = 5; # number of BZs in each dimension / 2

BZ_centers = reshape(collect(Base.product(G.*(-NBZ:NBZ), G.*(-NBZ:NBZ))), :)
BZ_kpoints = reshape(collect(Base.product(range(-G/2,G/2,Nk), range(-G/2,G/2,Nk))), :)

ϵ = Float32(G/1000) # numerical error tolerance

function get_Hk(
    BZ_kpoint::Tuple{Float32,Float32},
    BZ_centers::Vector{Tuple{Float32,Float32}},
    m::Float32,
    ϵ::Float32,
)
    n = length(BZ_centers)
    Hk = zeros(Float32, n, n)

    # diagonal
    pref = ħ^2 / e
    for i in 1:n
        kx = BZ_kpoint[1] + BZ_centers[i][1]
        ky = BZ_kpoint[2] + BZ_centers[i][2]
        Hk[i,i] = pref * (kx^2 + ky^2) / (2m)
    end

    # off-diagonal
    for i in 1:n, j in i+1:n
        if abs(norm(BZ_centers[i] .- BZ_centers[j]) - G) < ϵ
            Hk[i,j] = U0/2
            Hk[j,i] = U0/2
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

eta = (maximum(energies_LB)-minimum(energies_LB))/100         # broadening for DOS [eV]
kd = kde(energies_LB, bandwidth = eta)

# PLOTTING AND OUTPUT ====================================================
plt = plot(
    kd.x,
    kd.density,
    linewidth = 2,
    xlabel = "Energy",
    ylabel = "Density of States",
    yticks = false,
    legend = false,
    color = :red,
    framestyle = :box,
    title = "ZFDOS for U0=$(U0)eV, a=$(a*1f9)nm",
    size = (600,400),
    margin = 5mm,
)
png_filename = joinpath(outfolder_plots, "ZFDOS_U0=$(U0)_a=$(a*1f9).png")
savefig(plt, png_filename)

# save data
data_filename = joinpath(outfolder_data, "ZFDOS_U0=$(U0)_a=$(a*1f9).npy")
npzwrite(data_filename, energies_LB)