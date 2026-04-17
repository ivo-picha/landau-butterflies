# plots Bloch functions from the lowest hofstadter band
# for a given p, q, U0, X and Y
# outputs plots 

using Plots
using LinearAlgebra
using ProgressMeter
using Measures

include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis
include(joinpath(dirname(@__DIR__),"funcs/states.jl"))
using .States                        # build a Hamiltonian matrix in a Landau level basis

output_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/bloch_out"
mkpath(output_folder)

p = 1
q = 2
U0 = 0.05f0

pi32 = Float32(π)

X = 0.5f0*pi32
Y = 0.5f0*pi32
LLmax = 50

a_nm = 5.0 # lattice constant in nm
a = Float32(a_nm*1f-9) # in m
phi = Float32(p/q)

H = Hamil.get_full_ham(phi, X, Y, U0, a, p, LLmax)
evals, evecs = eigen(H)

Ngrid = 128

xrange = range(0f0, Float32(2q*a), length = Ngrid*q)
yrange = range(0f0, Float32(2a), length = Ngrid)

bloch_array = Array{ComplexF32}(undef, Ngrid, Ngrid*q, q)
for m in 1:q
    println("Calculating Bloch function for band $m...")
    @showprogress for (i,x) in enumerate(xrange)
        for (j,y) in enumerate(yrange)
            bloch_array[j,i,m] = States.get_eigenstate_from_evec(evecs[:,m], x, y, X, Y, p, q, a, LLmax)
        end
    end
    bloch_array[:,:,m] ./= maximum(abs.(bloch_array[:,:,m])) # normalize the bloch function for better plotting

    # plot the bloch function
    pltabs = Plots.plot(xlabel="x/a", ylabel="y/a", title="|ψ|", framestyle=:box, size =  (q*500,500))
    pltarg = Plots.plot(xlabel="x/a", ylabel="y/a", title="arg(ψ)", framestyle=:box, size =  (q*500,500))
    heatmap!(pltabs, xrange./a, yrange./a, abs.(bloch_array[:,:,m]), color=:viridis, aspect_ratio=1)
    heatmap!(pltarg, xrange./a, yrange./a, angle.(bloch_array[:,:,m]), color=:hsv, aspect_ratio=1)
    pltcombi = plot(pltabs, pltarg, layout = (2,1), size = (q*500,1000), title = "Band $m; X/π=$(X/π), Y/π=$(Y/π)", framestyle=:box)
    savefig(pltcombi, joinpath(output_folder, "bloch_p-$p-q-$q-U0-$U0-X-$X-Y-$Y-band$(m).png"))
end


