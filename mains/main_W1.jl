# calculates the spectrum of the lowest band at a given flux p/(q=1)
# Wannierizes the problem and outputs t and μ as a function of U0

using LinearAlgebra
using Plots
using ProgressMeter
using Base.Threads
using Measures

include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis
include(joinpath(dirname(@__DIR__),"funcs/states.jl"))
using .States                        # build a Hamiltonian matrix in a Landau level basis

#output folder
output_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/wannier_states/"
mkpath(output_folder)

p = 1
q = 1 # don't change
U0_list = collect(range(0.0151, 0.05, 32)) # in eV
a_nm = 5.0 # lattice constant in nm
NXY = 129 # number of k-points in each direction
LLmax = 25

a = Float32(a_nm*1f-9) # in m
phi = Float32(p/q)

# XY lists
X_list = collect(range(0f0, Float32(2π*q), length = NXY+1))[1:end-1]
Y_list = collect(range(0f0, Float32(2π), length = NXY+1))[1:end-1]
XY_list = reshape(collect(Base.product(X_list, Y_list)),:)

#make function that does for a given U later, now finish for a single U and test
U0 = -0.04f0

states = Tuple{Float32, Vector{ComplexF32}, Float32, Float32}[] # (energy, eigenvec, X, Y)
println("Calculating eigenstates for U0 = $U0 at ϕ = $p/$q...")
@showprogress @threads for (X,Y) in XY_list
    H = Hamil.get_full_ham(phi, X, Y, U0, a, p, LLmax)
    evals, evecs = eigen(H)
    for i in 1:q
        push!(states, (evals[i], evecs[:,i], X, Y))
    end
end

# sort states by energy
states = sort(states, by = first)

# compute and plot density in one magnetic unit cell for a given state
Ngrid = 128
xplotrange = range(-2.5*a*q, 2.5*a*q, Ngrid)
yplotrange = range(-2.5*a, 2.5*a, Ngrid)
dx = step(xplotrange)
dy = step(yplotrange)
xplotrange = Float32.(collect(xplotrange))
yplotrange = Float32.(collect(yplotrange))
# xyplotlist = reshape(collect(Base.product(xplotrange,yplotrange)),:)

# array that stores all eigenstates; wf_array[xidx,yidx,Xidx,Yidx] gives wf at (x,y) for state at (X,Y)
wf_array = Array{ComplexF32}(undef, Ngrid, Ngrid, NXY, NXY)

# fill wf_array
println("Calculating wavefunctions on grid...")
@showprogress @threads for state in states
    ens, evecs, Xs, Ys = state
    Xidx = findfirst(isequal(Xs), X_list)
    Yidx = findfirst(isequal(Ys), Y_list)
    for (i,x) in enumerate(xplotrange)
        for (j,y) in enumerate(yplotrange)
            wf_array[i,j,Xidx,Yidx] = States.get_eigenstate(x,y,state,p,q,a,LLmax)
        end
    end
end


# normalize densities
for (i,Xi) in enumerate(X_list)
    for (j,Yj) in enumerate(Y_list)
        wf_abs2 = abs2.(wf_array[:,:,i,j])
        norm = sum(wf_abs2) * dx * dy
        wf_array[:,:,i,j] .= wf_array[:,:,i,j] ./ sqrt(norm)
    end
end

# get trial wavefunction on a grid
σ = 0.3f0*a
g_array = Array{ComplexF32}(undef, Ngrid, Ngrid)
println("Calculating trial wavefunction on grid...")
for (i,x) in enumerate(xplotrange)
    for (j,y) in enumerate(yplotrange)
        g_array[i,j] = States.gaussian(x,y,0f0*a,0f0*a,σ)
    end
end

# get Loewdin Bloch functions
loewdin_array = Array{ComplexF32}(undef, Ngrid, Ngrid, NXY, NXY)
println("Calculating overlaps...")
@showprogress for (i,Xi) in enumerate(X_list)
    for (j,Yj) in enumerate(Y_list)
        wf = wf_array[:,:,i,j]
        overlap = sum(conj.(wf) .* g_array) * dx * dy
        phase = overlap / abs(overlap)
        loewdin_array[:,:,i,j] .= wf_array[:,:,i,j] .* phase
    end
end

# get wannier function vector
println("Calculating Wannier functions...")
Rx_list = Float32.(collect(-2q:1:2q))
Ry_list = Float32.(collect(-2:1:2))
Rxy_list = reshape(collect(Base.product(Rx_list, Ry_list)),:)
wannier_vector = Matrix{ComplexF32}[]
hwannier_vector = Matrix{ComplexF32}[]
@showprogress for (Rx, Ry) in Rxy_list
    wannier_array = States.get_wannier_array(q, -Rx, Ry, loewdin_array, X_list, Y_list)
    push!(wannier_vector, wannier_array)
    hw = States.H_on_wannier(wannier_array, xplotrange, yplotrange, phi, U0, a)
    push!(hwannier_vector, hw)
end

# get t and μ from wannier functions and save them
println("Calculating μ...")
mu_list = similar(Rxy_list, ComplexF32)
for i in eachindex(Rxy_list)
    # μ(r) = ⟨wᵣ|H|wᵣ⟩ / ⟨wᵣ|wᵣ⟩
    mu = sum(conj.(wannier_vector[i]) .* hwannier_vector[i])
    mu /= sum(abs2.(wannier_vector[i]))
    mu_list[i] = mu
    println("μ($(Rxy_list[i])) = $mu")
end

println("Calculating t...")
tx_array = zeros(ComplexF32, length(Rx_list)-1, length(Ry_list)) # hopping along x
ty_array = zeros(ComplexF32, length(Rx_list), length(Ry_list)-1) # along y
for (i,(Rxi,Ryi)) in enumerate(Rxy_list)
    for (j,(Rxj,Ryj)) in enumerate(Rxy_list)
        # only NN and positive direction of hopping
        if (Rxj,Ryj).-(Rxi,Ryi) == (1,0) # along x
            n = findfirst(x -> x == Rxi, Rx_list)
            m = findfirst(x -> x == Ryi, Ry_list)
            tx = sum(conj.(wannier_vector[i]) .* hwannier_vector[j])
            tx /= sum(conj.(wannier_vector[i]) .* wannier_vector[j])
            tx_array[n,m] = tx
        elseif (Rxj,Ryj).-(Rxi,Ryi) == (0,1) # along y
            n = findfirst(x -> x == Rxi, Rx_list)
            m = findfirst(x -> x == Ryi, Ry_list)
            ty = sum(conj.(wannier_vector[i]) .* hwannier_vector[j])
            ty /= sum(conj.(wannier_vector[i]) .* wannier_vector[j])
            ty_array[n,m] = ty
        end
    end
end

## out
println("\n|tₓ|:")
show(stdout, "text/plain", round.(abs.(tx_array); digits = 3))
println("\narg(tₓ):")
show(stdout, "text/plain", round.(angle.(tx_array); digits = 3))
println("\n|ty|:")
show(stdout, "text/plain", round.(abs.(ty_array); digits = 3))
println("\narg(ty):")
show(stdout, "text/plain", round.(angle.(ty_array); digits = 3))