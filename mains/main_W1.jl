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

p = 1
q = 1 # don't change
U0_list = collect(range(0.0151, 0.05, 32)) # in eV
a_nm = 5.0 # lattice constant in nm
NXY = 100 # number of k-points in each direction
LLmax = 20

a = Float32(a_nm*1f-9) # in m
phi = Float32(p/q)

# XY lists
X_list = collect(range(0f0, Float32(2π*q), length = NXY+1))[1:end-1]
Y_list = collect(range(0f0, Float32(2π), length = NXY+1))[1:end-1]
XY_list = reshape(collect(Base.product(X_list, Y_list)),:)

#make function that does for a given U later, now finish for a single U and test
U0 = -0.05f0

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
Ngrid = 64
xplotrange = range(-1.5*a*q,1.5*a*q,Ngrid)
yplotrange = range(-1.5*a*q,1.5*a,Ngrid)
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

# get wannier function array 
wannier_array = States.get_wannier_array(q, 0f0, 0f0, loewdin_array, X_list, Y_list)

# visualize wannier function
wf_real = real.(wannier_array)
wf_imag = imag.(wannier_array)
wf_abs2 = abs2.(wannier_array)
wf_arg = angle.(wannier_array)
plt_r = heatmap(xplotrange./a, yplotrange./a, wf_real', xlabel="x/a", ylabel="y/a", title="Re(ψ)", aspect_ratio=1, size=(500,500), legend = false);
plt_i = heatmap(xplotrange./a, yplotrange./a, wf_imag', xlabel="x/a", ylabel="y/a", title="Im(ψ)", aspect_ratio=1, size=(500,500), legend = false);
plt_a = heatmap(xplotrange./a, yplotrange./a, sqrt.(wf_abs2)', xlabel="x/a", ylabel="y/a", title="|ψ|", aspect_ratio=1, size=(500,500), legend = false);
plt_th = heatmap(xplotrange./a, yplotrange./a, reshape(wf_arg', (Ngrid, Ngrid))', xlabel="x/a", ylabel="y/a", title="arg(ψ)", aspect_ratio=1, size=(500,500), legend = false);
plt_combo = plot(plt_r, plt_i, plt_a, plt_th, layout = grid(2, 2, hgap = 2mm, vgap = 2mm),
    aspect_ratio = 1,
    size = (800,800),
    margin = 1mm, plot_title="Wannier; ϕ=$(p)/$(q), U0=$(U0) eV")
