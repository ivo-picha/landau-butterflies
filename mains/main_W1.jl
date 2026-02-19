# calculates the spectrum of the lowest band at a given flux p/(q=1)
# Wannierizes the problem and outputs t and μ as a function of U0

using LinearAlgebra
using Plots
using ProgressMeter
using Base.Threads
using Measures
using DelimitedFiles

include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis
include(joinpath(dirname(@__DIR__),"funcs/states.jl"))
using .States                        # build a Hamiltonian matrix in a Landau level basis

p = 2
q = 1 # don't change
a_nm = 5.0 # lattice constant in nm
NXY = 256 # number of k-points in each direction
LLmax = 20

a = Float32(a_nm*1f-9) # in m
phi = Float32(p/q)

# XY lists
X_list = collect(range(0f0, Float32(2π*q), length = NXY+1))[1:end-1]
Y_list = collect(range(0f0, Float32(2π), length = NXY+1))[1:end-1]
XY_list = reshape(collect(Base.product(X_list, Y_list)),:)

#make function that does for a given U later, now finish for a single U and test
U0 = -0.05f0

#output folder
output_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/wannier_out/ph$p-$q-U$U0-a$a_nm-LL$LLmax-NXY$NXY"
mkpath(output_folder)

states = Tuple{Float32, Vector{ComplexF32}, Float32, Float32}[] # (energy, eigenvec, X, Y)
println("Calculating eigenstates for U0 = $U0 at ϕ = $p/$q...")
@showprogress for (X,Y) in XY_list
    H = Hamil.get_full_ham(phi, X, Y, U0, a, p, LLmax)
    evals, evecs = eigen(H)
    for i in 1:q
        push!(states, (evals[i], evecs[:,i], X, Y))
    end
end

# sort states by energy
states = sort(states, by = first)

# extract energies and (X,Y) for FT; "states" will be deleted later
energies = [];
XYs = [];
for state in states
    push!(energies, state[1])
    push!(XYs, (state[3],state[4]))
end

# compute and plot density in one magnetic unit cell for a given state
Ngrid = 128
xplotrange = range(-1.5*a*q, 1.5*a*q, Ngrid)
yplotrange = range(-1.5*a, 1.5*a, Ngrid)
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

states = nothing # free memory

# normalize densities
for (i,Xi) in enumerate(X_list)
    for (j,Yj) in enumerate(Y_list)
        wf_abs2 = abs2.(wf_array[:,:,i,j])
        norm = sum(wf_abs2) * dx * dy
        wf_array[:,:,i,j] .= wf_array[:,:,i,j] ./ sqrt(norm)
    end
end

states = nothing # free memory


# get trial wavefunction on a grid
g_array = Array{ComplexF32}(undef, Ngrid, Ngrid)
println("Calculating trial wavefunction on grid...")
for (i,x) in enumerate(xplotrange)
    for (j,y) in enumerate(yplotrange)
        g_array[i,j] = States.gaussian_LG(x,y,0f0*a,0f0*a,phi,U0,a) # centered at (0,0); works with other UC centres
    end
end

# get Loewdin Bloch functions
println("Calculating overlaps (Loewdin functions)...")
@showprogress for (i,Xi) in enumerate(X_list)
    for (j,Yj) in enumerate(Y_list)
        wf = wf_array[:,:,i,j]
        overlap = sum(conj.(wf) .* g_array) * dx * dy
        phase = overlap / abs(overlap)
        wf_array[:,:,i,j] .*= phase
    end
end

# get wannier function vectorLoewd
println("Calculating Wannier functions...")
Rx_list = Float32.(collect(-1q:1:1q))
Ry_list = Float32.(collect(-1:1:1))
Rxy_list = reshape(collect(Base.product(Rx_list, Ry_list)),:)
wannier_vector = Matrix{ComplexF32}[];
# hwannier_vector = Matrix{ComplexF32}[]
@showprogress for (j,(Rx, Ry)) in enumerate(Rxy_list)
    wannier_array = States.get_wannier_array(q, Rx, Ry, wf_array, X_list, Y_list)
    push!(wannier_vector, wannier_array ./ sqrt(dx * dy)) # cleaning up normalization
end

wf_array = nothing # free memory

# save wannier_vector, hwannier_vector, and plotranges as npz
# npzwrite(joinpath(output_folder, "wannier_output_ph$p-$q-U$U0-NXY$NXY-Ng$Ngrid.npz"), Dict(
#     "wannier_vector" => wannier_vector,
#     "hwannier_vector" => hwannier_vector,
#     "xplotrange" => xplotrange,
#     "yplotrange" => yplotrange,
#     "Rx_list" => Rx_list,
#     "Ry_list" => Ry_list
# ))




## visualize all wannier functions on one plot

wannier_superarray = zeros(ComplexF32, Ngrid, Ngrid)
for (i,(Rx,Ry)) in enumerate(Rxy_list)
    wannier_superarray .+= wannier_vector[i]
end
wannier_superarray ./= maximum(abs.(wannier_superarray))
pltsize = (476,500)
plt_real = heatmap(xplotrange./a, yplotrange./a, real.(wannier_superarray)', xlabel="x/a", ylabel="y/a", title="Re(sum(wᵣ))", aspect_ratio=1, size=pltsize)
plt_imag = heatmap(xplotrange./a, yplotrange./a, imag.(wannier_superarray)', xlabel="x/a", ylabel="y/a", title="Im(sum(wᵣ))", aspect_ratio=1, size=pltsize)
plt_abs = heatmap(xplotrange./a, yplotrange./a, abs.(wannier_superarray)', xlabel="x/a", ylabel="y/a", title="|sum(wᵣ)|", aspect_ratio=1, size=pltsize)
plt_arg = heatmap(xplotrange./a, yplotrange./a, angle.(wannier_superarray)', xlabel="x/a", ylabel="y/a", title="arg(sum(wᵣ))", aspect_ratio=1, size=pltsize)

hline!(plt_real, [-1.5,-0.5,0.5,1.5], color=:white, linestyle=:dash, label = "")
vline!(plt_real, [-1.5,-0.5,0.5,1.5], color=:white, linestyle=:dash, label = "")
hline!(plt_imag, [-1.5,-0.5,0.5,1.5], color=:white, linestyle=:dash, label = "")
vline!(plt_imag, [-1.5,-0.5,0.5,1.5], color=:white, linestyle=:dash, label = "")
hline!(plt_abs, [-1.5,-0.5,0.5,1.5], color=:white, linestyle=:dash, label = "")
vline!(plt_abs, [-1.5,-0.5,0.5,1.5], color=:white, linestyle=:dash, label = "")
hline!(plt_arg, [-1.5,-0.5,0.5,1.5], color=:white, linestyle=:dash, label = "")
vline!(plt_arg, [-1.5,-0.5,0.5,1.5], color=:white, linestyle=:dash, label = "")

savefig(plt_real, joinpath(output_folder, "wannier_superarray_real_ph$p-$q-U$U0.png"))
savefig(plt_imag, joinpath(output_folder, "wannier_superarray_imag_ph$p-$q-U$U0.png"))
savefig(plt_abs, joinpath(output_folder, "wannier_superarray_abs_ph$p-$q-U$U0.png"))
savefig(plt_arg, joinpath(output_folder, "wannier_superarray_arg_ph$p-$q-U$U0.png"))
println("Saved superarray plots to $output_folder")







## ANALYSIS OF WANNIER FUNCTIONS

# get overlap of a specified H.wannier function with all wannier functions; make a matrix

function get_wHw_matrix(wR::Tuple{Real,Real})

    wannierR = wannier_vector[findfirst(x->x==wR, Rxy_list)]
    hwannierR = States.H_on_wannier(wannierR, xplotrange, yplotrange, phi, U0, a)

    whw_matrix = zeros(ComplexF32, length(Rx_list), length(Ry_list))
    for (i,(Rxi,Ryi)) in enumerate(Rxy_list)
            whw = sum(conj.(wannier_vector[i]) .* hwannierR)

            n = findfirst(x -> x == Rxi, Rx_list)
            m = findfirst(x -> x == Ryi, Ry_list)
            whw_matrix[n,m] = whw 
    end
    #whw_matrix ./= maximum(abs.(whw_matrix))
    #show(stdout, "text/plain", whw_matrix)
    # println("Overlap matrix of wᵣ with H.w₁,₀; normalized:")
    # show(stdout, "text/plain", round.(whw_matrix; digits = 5))


    # energy FT

    t_ft_matrix = zeros(ComplexF32, length(Rx_list), length(Ry_list))
    for (i,(Rxi,Ryi)) in enumerate(Rxy_list)
        tR = 0f0 + im*0f0
        for (i,(X,Y)) in enumerate(XYs)
            tR += energies[i] * exp(-im*(X*(Ryi-wR[2]) - Y*(Rxi-wR[1]))/q)
        end
        n = findfirst(x -> x == Rxi, Rx_list)
        m = findfirst(x -> x == Ryi, Ry_list)
        t_ft_matrix[n,m] = tR * q / NXY^2
    end
    #t_ft_matrix ./= maximum(abs.(t_ft_matrix))
    # println("tR from energy FT; normalized:")
    # show(stdout, "text/plain", t_ft_matrix; digits = 5)

    return whw_matrix, t_ft_matrix
end

whw_matrix00, t_ft_matrix00 = get_wHw_matrix((0f0,0f0))
whw_matrix10, t_ft_matrix10 = get_wHw_matrix((1f0,0f0))
whw_matrix01, t_ft_matrix01 = get_wHw_matrix((0f0,1f0)) 

whw_tx_rescaled = whw_matrix00[3,2] * t_ft_matrix00[2,2] / whw_matrix00[2,2]
whw_ty_rescaled = whw_matrix00[2,3] * t_ft_matrix00[2,2] / whw_matrix00[2,2]
# save matrices as one text file in folder
whw_matrices_path = joinpath(output_folder, "whw_and_ft_matrices.txt")

matrices = transpose.([whw_matrix00, t_ft_matrix00, whw_matrix10, t_ft_matrix10, whw_matrix01, t_ft_matrix01])
descriptions = [
    "^\n|\ny\nx-->\n\n Center at (0,0)\nFrom overlap <wᵣ|H|w₀,₀>\n",
    "Rescaled tₓ: $whw_tx_rescaled\nRescaled t_y: $whw_ty_rescaled\nFrom energy FT:\n",
    "\n\nCenter at (1,0)\nFrom overlap <wᵣ|H|w₁,₀>t_ft_matrix = zeros(ComplexF32, length(Rx_list), length(Ry_list))
    for (i,(Rxi,Ryi)) in enumerate(Rxy_list)
        tR = 0f0 + im*0f0
        for (i,(X,Y)) in enumerate(XYs)
            tR += energies[i] * exp(-im*(X*(Ryi-wR[2]) - Y*(Rxi-wR[1]))/q)
        end
        n = findfirst(x -> x == Rxi, Rx_list)
        m = findfirst(x -> x == Ryi, Ry_list)
        t_ft_matrix[n,m] = tR * q / NXY^2
    end\n",
    "From energy FT:\n",
    "\n\nCenter at (0,1)\nFrom overlap <wᵣ|H|w₀,₁>\n",
    "From energy FT:\n",
]

function write_matrix_aligned(io, mat)
    rows, cols = size(mat)
    # Determine column width based on the largest number in the matrix
    max_width = maximum(length.(string.(mat))) + 1  # +1 for spacing

    # Optional: add column headers
    header = "     " * join([rpad("C$j", max_width) for j in 1:cols])
    write(io, header * "\n")

    for i in 1:rows
        line = rpad("R$i", 5) * join([rpad(string(mat[i,j]), max_width) for j in 1:cols])
        write(io, line * "\n")
    end
end
open(whw_matrices_path, "w") do io
    for (desc, mat) in zip(descriptions, matrices)
        write(io, desc * "\n")         # write the description
        writedlm(io, mat)              # write the matrix
        write(io, "\n")                # extra line for spacing
    end
end
println("Saved whw and t_ft matrices to $whw_matrices_path")
# Save all matrices to one file
open(whw_matrices_path, "w") do io
    for (desc, mat) in zip(descriptions, matrices)
        write(io, desc * "\n")
        write_matrix_aligned(io, mat)
        write(io, "\n")
    end
end
println("Saved whw and tR from FT matrices to $whw_matrices_path")

## 
