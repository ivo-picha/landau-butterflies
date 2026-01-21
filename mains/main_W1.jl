# calculates the spectrum of the lowest band at a given flux p/(q=1)
# Wannierizes the problem and outputs t and μ as a function of U0

using LinearAlgebra
using Plots
using ProgressMeter
using Base.Threads

include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis
include(joinpath(dirname(@__DIR__),"funcs/states.jl"))
using .States                        # build a Hamiltonian matrix in a Landau level basis

p = 1
q = 1 # don't change
U0_list = collect(range(0.0151, 0.05, 32)) # in eV
a_nm = 5.0 # lattice constant in nm
NXY = 65 # number of k-points in each direction
LLmax = 30

a = Float32(a_nm*1f-9) # in m
phi = Float32(p/q)

# XY lists
X_list = collect(range(0f0, Float32(2π*q), length = NXY+1))[1:end-1]
Y_list = collect(range(0f0, Float32(2π), length = NXY+1))[1:end-1]
XY_list = reshape(collect(Base.product(X_list, Y_list)),:)

#make function that does for a given U later, now finish for a single U and test
U0 = 0.04f0

states = Tuple{Float32, Vector{ComplexF32}, Float32, Float32}[] # (energy, eigenvec, X, Y)

@showprogress for (X,Y) in XY_list
    H = Hamil.get_full_ham(phi, X, Y, U0, a, p, LLmax)
    evals, evecs = eigen(H)
    for i in 1:q
        push!(states, (evals[i], evecs[:,i], X, Y))
    end
end

# sort states by energy
states = sort(states, by = first)

Ngrid = 32
xplotrange = collect(range(0,5*a*q,Ngrid))
yplotrange = collect(range(0,5*a,Ngrid))
xyplotlist = reshape(collect(Base.product(xplotrange,yplotrange)),:)

dens = [];
ss = 1
@showprogress for (x,y) in xyplotlist
    push!(dens, abs(States.get_eigenstate(x,y,states[ss],p,q,a,LLmax)))
end

dens_grid = reshape(dens, (Ngrid, Ngrid))
heatmap(xplotrange./a, yplotrange./a, dens_grid', xlabel="x", ylabel="y", title="Density of state #$(ss) at U0=$(U0) eV", aspect_ratio=1)




