# spectrum and states for 2d electron gas in cos potential
using LinearAlgebra
using ProgressMeter
using Plots

plot_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/densities_zerofield"
#plot_folder = "/users/ivoga/lh-short/plts"

#args = ARGS
args = ["[10, 4]"]

# parse arguments
args1 = replace(args[1], "[" => "", "]" => "")
args2 = split(args1, ",")
args_vec = parse.(Float64, args2)

nu = Float64(args_vec[2]); # filling factor
a = 1.; # lattice constant
U0 = Float64(args_vec[1]); # potential strength
m = 1.; # electron/particle mass

kF = sqrt(4π*nu)/a # Fermi wavevector (spinless electrons in 2D)
G = 2π/a # recip scat vec

Nk = 30; # number of momentum states in 1D of the BZ

k1_list = range(-kF, kF, step = G/Nk)
k_grid = reshape(collect(Base.product(k1_list, k1_list)), :) # grid in momentum space

filter!(k -> ((k[1])^2 + (k[2])^2) <= kF^2, k_grid) # only take states below Fermi; also defines a basis

# building the hamiltonian
println("Building a $(length(k_grid)) x $(length(k_grid)) Hamiltonian matrix.")
E_kin(k) = ((k[1])^2 + (k[2])^2)/(2*m) # kinetic energy; ħ=1
H = diagm(E_kin.(k_grid)) # diagonal part of the matrix

# add scattering from potential
eps = (G/Nk)/1000 # numerical error tolerance
@showprogress for n = eachindex(k_grid)
    for m = eachindex(k_grid)
        if ((abs(k_grid[n][1] - (k_grid[m][1] + G)) < eps) &&  (abs(k_grid[n][2] - k_grid[m][2]) < eps))
            H[n,m] = U0
            H[m,n] = U0
        elseif ((abs(k_grid[n][2] - (k_grid[m][2] + G)) < eps) &&  (abs(k_grid[n][1] - k_grid[m][1]) < eps))
            H[n,m] = U0
            H[m,n] = U0
        end
    end
end


if ishermitian(H)
    H = Hermitian(H)
else
    println("Error: Hamiltonian is not Hermitian.")
    exit(1)
end

# diagonalize
evals, evecs = eigen(H)
evecsR = round.(evecs; digits = 10) # ease further calculations

# wavefunction functions
function wf(x::Float64, y::Float64, kx::Float64, ky::Float64)::ComplexF64
    return exp(im*(x*kx + y*ky))    
end

function get_density(x::Float64, y::Float64, k_grid::Vector{Tuple{Float64, Float64}}, evecs::Matrix{Float64})::Float64
    tot = 0.0;
    for i = 1:Int(sqrt(length(evecsR)))
        vec = evecs[:,i]
        psi = sum([ vec[n]*wf(x,y,k_grid[n][1],k_grid[n][2]) for n = eachindex(k_grid)])
        tot = tot + psi*conj(psi)
    end
    return real(tot)/(Int(sqrt(length(evecsR))))
end


# create a grid in real space
Nr = 100; # number of steps in total
Nuc = 3; # number of unit cells plotted in each dim
r_range = range(0.,a*Nuc, Nr)
r_grid = reshape(collect(Base.product(r_range, r_range)), :)

density_list = [];
println("Calculating real space denity with $(length(r_grid)) points.")
@showprogress for xy in r_grid
    push!(density_list, get_density(xy[1],xy[2],k_grid,evecs))
end
density_mat = transpose(reshape(density_list, Nr, Nr))

# plot
p1 = heatmap(collect(r_range), collect(r_range), density_mat)
title!(p1, "U₀=$U0, a=$a, ν=$nu")

# save plot
plot_name = "ZF_U0$U0-a$a-nu$nu.png"
savefig(p1, joinpath(plot_folder,plot_name))