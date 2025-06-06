# spectrum and states for 2d electron gas in cos potential
using LinearAlgebra
using ProgressMeter
using Plots
using Statistics
using NPZ

using Base.Threads
nt = nthreads()

plot_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/densities_zerofield"
#plot_folder = "/users/ivoga/lh-short/plts"

const ħ = 6.62607015e-34/(2π);  # Planck constant [J s]
const e = 1.602176634e-19;      # elementary charge [C]
const m_e = 9.1093837139e-31;   # electron mass [kg];

#args = ARGS
args = ["[0.01, 1.2]"] # [U0 (eV), filling]

# parse arguments
args1 = replace(args[1], "[" => "", "]" => "")
args2 = split(args1, ",")
args_vec = parse.(Float64, args2)

nu = Float64(args_vec[2]); # filling factor
a = 5e-9; # lattice constant
U0 = Float64(args_vec[1]); # potential strength in eV
m = m_e; # electron/particle mass

G = 2π/a # recip scat vec
EF = (ħ^2 /e)* 2π*(nu/a^2)/(2*m) 

Nk = 33; # sqrt of number of momentum states in BZ
NBZ = 10; # number of BZs in each dimension / 2

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
                Hk[i,j] = U0
            end
        end
    end
    return Hermitian(Hk)
end

println("Calculating spectrum and eigenstates...")

# store state information like (energy, (kx, ky) ∈ BZ , eigenvector)
#states = Tuple{Float64, Tuple{Float64, Float64}, Vector{ComplexF64}}[];
states_buffers = [Vector{Tuple{Float64, Tuple{Float64, Float64}, Vector{ComplexF64}}}() for _ in 1:nt]
@showprogress for kpoint in BZ_kpoints
    Hk = get_Hk(kpoint, BZ_centers, m, ϵ)
    evals, evecs = eigen(Hk)
    @threads for i = eachindex(evals)
        tid = threadid()
        push!(states_buffers[tid], (evals[i], kpoint, evecs[:,i]))
    end
end

states = reduce(vcat, states_buffers)

# discard states above Fermi level
states_sorted = sort(states, by = first)
states_filtered = filter(s -> first(s) <= EF, states_sorted)


# wavefunction functions
function wf(x::Float64, y::Float64, kx::Float64, ky::Float64)::ComplexF64
    return exp(im*(x*kx + y*ky))    
end

function get_density(x::Float64, y::Float64, BZ_centers::Vector{Tuple{Float64, Float64}}, states_filtered::Vector{Tuple{Float64, Tuple{Float64, Float64}, Vector{ComplexF64}}})::Float64
    tot = 0.0;
    
    for state in states_filtered
        psi = sum([ state[3][i] * wf(x,y,(state[2][1] + BZ_centers[i][1]), (state[2][2] + BZ_centers[i][2])) for i = eachindex(BZ_centers)])
        tot += psi*conj(psi)
    end

    return real(tot)
end


# create a grid in real space
Nr = 32; # number of steps in total
Nuc = 1; # number of unit cells plotted in each dim
r_range = range(0.,a*Nuc, Nr)
r_grid = reshape(collect(Base.product(r_range, r_range)), :)

density_list = Float64[];
println("Calculating real space denity with $(length(r_grid)) points...")
@showprogress for xy in r_grid
    push!(density_list, get_density(xy[1],xy[2],BZ_centers,states_filtered))
end

norm_factor = nu/mean(density_list) #/a^2
density_mat = transpose(reshape(norm_factor.*density_list, Nr, Nr))

# plot
p1 = heatmap(collect(r_range), collect(r_range), density_mat, clim = (0,2*nu))
title!(p1, "U₀=$U0, a=$a, ν=$nu")

# save data so it can be accessed later; npz format, readable by python as well
npzwrite(joinpath(plot_folder, "ZF_U0$U0-a$a-nu$nu.npz"),
        Dict("x" => collect(r_range), "y" => collect(r_range), "z" => density_mat)
)

# save plot
plot_name = "ZF_U0$U0-a$a-nu$nu.png"
savefig(p1, joinpath(plot_folder,plot_name))
