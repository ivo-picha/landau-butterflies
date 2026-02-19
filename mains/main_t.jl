# calculates the hopping amplitudes for a given U₀ and ϕ=p/q
# FT approach, no Wannier construction, just overlap matrices w/ trial WFs, setting a gauge

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

p = 1
q = 2 
U0 = -0.05f0 # keep negative (energy minimum at origin)
a_nm = 5.0 # lattice constant in nm
NXY = 256 # number of k-points in each direction
LLmax = q*25

a = Float32(a_nm*1f-9) # in m
phi = Float32(p/q)

# XY lists
X_list = collect(range(0f0, Float32(2π*q), length = NXY+1))[1:end-1]
Y_list = collect(range(0f0, Float32(2π), length = NXY+1))[1:end-1]
XY_list = reshape(collect(Base.product(X_list, Y_list)),:)

#output folder
output_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/wannier_out/ph$p-$q-U$U0-a$a_nm-LL$LLmax-NXY$NXY"
# cluster path
#output_folder = "/users/ivoga/lh/out/wannier_out/ph$p-$q-U$U0-a$a_nm-LL$LLmax-NXY$NXY"
mkpath(output_folder)


# real space grid
Ngrid = 100
x_grid = Float32.(collect(range(-a*(q+0.5), a*(q+0.5), length = Ngrid*q)))
y_grid = Float32.(collect(range(-1.5*a, 1.5*a, length = Ngrid)))

# trial wavefunctions
println("Generating trial wavefunctions...")
g_array = Array{ComplexF32}(undef, Ngrid*q, Ngrid, q)
normg = 0f0
for m in 1:q
    for (i,x) in enumerate(x_grid)
        for (j,y) in enumerate(y_grid)
            g_xy = States.gaussian_LG(x,y,Float32(m-1)*a,0f0*a,phi,U0,a)
            g_array[i,j,m] = g_xy
            normg += abs2(g_xy)
        end
    end
    global normg
    normg *= (x_grid[2] - x_grid[1])*(y_grid[2] - y_grid[1])
    normg = sqrt(normg)
    g_array[:,:,m] ./= normg
end

# loop over XY
C_array = Array{ComplexF32}(undef, q, q, NXY, NXY)
println("Calculating eigenstates for U0 = $(abs(U0)) at ϕ = $p/$q...")

@showprogress @threads for i = 1:NXY
    X = X_list[i]
    for j = 1:NXY
        Y = Y_list[j]

        states = Tuple{Float32, Vector{ComplexF32}}[] # (energy, eigenvec, X, Y)
        energies = Float32[];

        H = Hamil.get_full_ham(phi, X, Y, U0, a, p, LLmax)
        evals, evecs = eigen(H)
        for i in 1:q
            push!(states, (evals[i], evecs[:,i]))
            push!(energies, evals[i])
        end

        # sort states by energy
        sort!(energies)
        sort!(states, by = first)

        # array to store wavefunctions
        wf_array = Array{ComplexF32}(undef, Ngrid*q, Ngrid, q)
        norm = 0f0
        for m in 1:q
            for (i,x) in enumerate(x_grid)
                for (j,y) in enumerate(y_grid)
                    wf_xy = States.get_eigenstate_XY(x,y,states[m],X,Y,p,q,a,LLmax)
                    wf_array[i,j,m] = wf_xy
                    norm += abs2(wf_xy)
                end
            end
            norm *= (x_grid[2] - x_grid[1])*(y_grid[2] - y_grid[1])
            norm = sqrt(norm)
            wf_array[:,:,m] ./= norm
        end

        # matrix Aₘₙ=⟨ψₘ|gₙ⟩
        Amat = Array{ComplexF32}(undef, q, q)
        for m in 1:q
            for n in 1:q
                Amn = sum(g_array[:,:,m] .* conj.(wf_array[:,:,n])) * (x_grid[2] - x_grid[1])*(y_grid[2] - y_grid[1])
                Amat[m,n] = Amn
            end
        end

        # overlap matrix S = A†.A
        Smat = Hermitian(Amat' * Amat)
        Smat_inv_sqrt = inv(sqrt(Smat))
        ASinv_sqrt = Amat * Smat_inv_sqrt
        Cmat = ASinv_sqrt' * diagm(energies) * ASinv_sqrt

        C_array[:,:,i,j] = Cmat
    end
end

# wR: Wannier centre; set to zero for FT approach
# nR: wannier state number (within unit cell); up to q
function get_hops(wR::Tuple{Real,Real},nR::Integer)
    Rxy_list = [wR .+ ((i-2)*1f0,(j-2)*1f0) for i in 1:3 for j in 1:3] # list of 8 NN coordinates = 4 NNs + 4 diagonals; + origin for chem potential
    # list of wannier state numbers at the coordinates above
    mR_list = [mod1(nR-1,q), mod1(nR,q), mod1(nR+1,q), mod1(nR-1,q), mod1(nR,q), mod1(nR+1,q), mod1(nR-1,q), mod1(nR,q), mod1(nR+1,q)]

    t_ft_list = ComplexF32[];
    for (n,(Rxn,Ryn)) in enumerate(Rxy_list)
        tR = 0f0 + im*0f0
        m = mR_list[n]
        for (i,X) in enumerate(X_list)
            for (j,Y) in enumerate(Y_list)
                tR += C_array[nR,m,i,j] * exp(-im*(X*(wR[2]-Ryn) - Y*(wR[1]-Rxn))/q)
            end
        end
        t_ft_list = push!(t_ft_list, tR / NXY^2)
    end

    absmat = reshape(abs.(t_ft_list), 3, 3)
    argmat = reshape(round.(mod.(angle.(-1.0 .* t_ft_list)/(2π) .+0.5 ,1.0) .- 0.5;digits = 5), 3, 3)
    return absmat, argmat
end

# function to write the matrix in aligned columns with row and column headers
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

wR_list = [(Float32(i), 0f0) for i = 0:(q-1)]
nR_list = collect(1:q)

# output
out_path = joinpath(output_folder, "whw_and_ft_matrices.txt")
open(out_path, "w") do io
    write(io, "Matrices for hopping amplitudes and chemical potentials on/to different sites.\n")
    write(io, "Each Wannier centre is at the centre of the matrix. Others are one lattice constant apart.\n")
    write(io, "x↓ y→\n")
    write(io, "================\n\n")
    for nR in nR_list
        wR = wR_list[nR]
        absmat, argmat = get_hops(wR, nR)
        write(io, "Wannier centre: $wR, Wannier state number: $nR\n")
        write(io, "Amplitudes |t|:\n")
        write_matrix_aligned(io, absmat)
        write(io, "Phases θ/2π:\n")
        write_matrix_aligned(io, argmat)
        write(io, "\n\n")
    end
end

println("Done! Hopping matrices written to $out_path.")