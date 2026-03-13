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

args = ARGS
if length(ARGS) != 4
    println("USAGE: main_t.jl p q U0 LLmax")
end

p = parse(Int, args[1])
q = parse(Int, args[2])
U0 = -parse(Float32, args[3]) # so that there is a wannier state at the origin
LLmax = parse(Int, args[4])


a_nm = 5.0 # lattice constant in nm
NXY = 20 # number of k-points in each direction; for larger p and q consider using this for every 2pi in X


a = Float32(a_nm*1f-9) # in m
phi = Float32(p/q)

# XY lists
X_list = collect(range(0f0, Float32(2π*q), length = NXY+1))[1:end-1]
Y_list = collect(range(0f0, Float32(2π), length = NXY+1))[1:end-1]
XY_list = reshape(collect(Base.product(X_list, Y_list)),:)

#output folder
output_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/wannier_out"
# cluster path
#output_folder = "/users/ivoga/lh/out/wannier_out"

mkpath(output_folder)

# real space grid
Ngrid = 20
x_grid = Float32.(collect(range(-a*(q+0.5), a*(q+0.5), length = Ngrid*q)))
y_grid = Float32.(collect(range(-1.5*a, 1.5*a, length = Ngrid)))

# trial wavefunctions
println("Generating trial wavefunctions...")
g_array = Array{ComplexF32}(undef, Ngrid*q, Ngrid, q)
normg = 0f0
for m in 1:q
    global normg
    for (i,x) in enumerate(x_grid)
        for (j,y) in enumerate(y_grid)
            g_xy = States.gaussian_LG(x,y,Float32(m-1)*a,0f0*a,phi,U0,a)
            g_array[i,j,m] = g_xy
            normg += abs2(g_xy)
        end
    end
    normg *= (x_grid[2] - x_grid[1])*(y_grid[2] - y_grid[1])
    normg = sqrt(normg)
    g_array[:,:,m] ./= normg
end

# loop over XY
C_array = Array{ComplexF32}(undef, q, q, NXY, NXY)
println("Calculating eigenstates for U₀ = $(abs(U0)) at ϕ = $p/$q...")

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
        # sort!(energies)
        sort!(states, by = first)

        # array to store wavefunctions
        wf_array = Array{ComplexF32}(undef, Ngrid*q, Ngrid, q)
        for m in 1:q
            norm = 0f0
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
                Amn = sum(g_array[:,:,n] .* conj.(wf_array[:,:,m])) * (x_grid[2] - x_grid[1])*(y_grid[2] - y_grid[1])
                Amat[m,n] = Amn
            end
        end

        #test
        # svds = svdvals(Amat)
        # if any(svds .< 1e-3)
        #     println("Warning: SVD values of A matrix are near zero at X=$X, Y=$Y")
        # end

        # overlap matrix S = A†.A
        Smat = Hermitian(Amat' * Amat)
        #Smat_inv_sqrt = inv(sqrt(Smat))
        #Umat = Amat * Smat_inv_sqrt
        #
        F = eigen(Hermitian(Smat))
        λ = real.(F.values)

        tol = 1e-12 * maximum(λ)
        if minimum(λ) < tol
            @warn "Overlap matrix nearly singular" minimum(λ) maximum(λ)
        end

        inv_sqrt_λ = [x > tol ? 1 / sqrt(x) : 0.0 for x in λ]
        Smat_inv_sqrt = F.vectors * Diagonal(inv_sqrt_λ) * F.vectors'
        Umat = Amat * Smat_inv_sqrt
        #
        Cmat = Umat' * diagm(energies) * Umat

        #test
        # Ucheck = norm(Umat'*Umat - I)
        # if Ucheck > 1e-4
        #     println("Warning: U matrix is not unitary at X=$X, Y=$Y, norm: $Ucheck")
        # end


        C_array[:,:,i,j] = Cmat
    end
end



# wR: Wannier centre; set to zero for FT approach
# nR: wannier state number (within unit cell); up to q
function get_hops(wR::Tuple{Real,Real},nR::Integer)
    Rxy_list = [wR .+ ((j-2)*1f0,(i-2)*1f0) for i in 1:3 for j in 1:3] # list of 8 NN coordinates = 4 NNs + 4 diagonals; + origin for chem potential
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
    absmat[2,2] = real(reshape(t_ft_list, 3, 3)[2,2]) # set chem potential to the actual value, not the abs value
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


# test
# tR = 0f0 + im*0f0
# n = 2
# wR = (0f0, 0f0)
# Rxn, Ryn = (0f0, 0f0)
# m = 1
# for (i,X) in enumerate(X_list)
#     for (j,Y) in enumerate(Y_list)
#         tR += C_array[n,m,i,j] * exp(-im*(X*(wR[2]-Ryn) - Y*(wR[1]-Rxn))/q)
#     end
# end
# tR /= NXY^2
# println("abs(tR) = ", abs(tR))
# println("angle(tR)/(2π) = ", angle(-tR)/(2π))






wR_list = [(Float32(i), 0f0) for i = 0:(q-1)]
nR_list = collect(1:q)

# output
out_path1 = joinpath(output_folder, "testhopping_amps_mats-ph$p-$q-U$U0-a$a_nm-LL$LLmax-NXY$NXY.txt")
open(out_path1, "w") do io
    write(io, "Matrices for hopping amplitudes and chemical potentials on/to different sites.\n")
    write(io, "Each Wannier centre is at the centre of the matrix. Others are one lattice constant apart.\n")
    write(io, "x↓ y→\n")
    write(io, "U₀=$(abs(U0)), ϕ=$p/$q, a=$(a_nm)nm, LLmax=$LLmax, NXY=$NXY, Ngrid=$Ngrid\n")
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

println("Done! Hopping matrices written to $out_path1.")


# write output in a wannier90 style: every hopping is included 
out_path2 = joinpath(output_folder, "all_hops-ph$p-$q-U$U0-a$a_nm-LL$LLmax-NXY$NXY.txt")
open(out_path2, "w") do io
    write(io, "Hopping amplitudes and chemical potentials (eV) for:\nU₀=$(abs(U0)), ϕ=$p/$q, a=$(a_nm)nm, LLmax=$LLmax, NXY=$NXY, Ngrid=$Ngrid.\n\n")
    write(io, "----------------------------\nChemical potentials\n----------------------------\n")
    write(io, "Rx   Ry  μ[eV]\n")
    #Calculate chemical potentials
    for m in 1:q
        mu_m = sum(C_array[m,m,:,:])
        write(io, "$(m-1)    0    $(real(mu_m))\n")
    end
    write(io, "\n\n----------------------------\nHopping amplitudes\n----------------------------\n")
    write(io, "From: Rx  Ry  n  To: Rx' Ry'  m          |t|[eV]      arg(t)/2π\n")
    Rto_list = [(ux+m,uy) for ux=-1:1 for uy=-2:2 for m=0:q-1]
    Rfrom_list = [(m,0) for m=0:q-1]
    for Rfrom in Rfrom_list
        n = Rfrom[1]
        for Rto in Rto_list
            m = mod(Rto[1],q)
            Rto == Rfrom && continue
            t = 0f0 + im*0f0
            for (i,X) in enumerate(X_list)
                for (j,Y) in enumerate(Y_list)
                    t += C_array[n+1,m+1,i,j] * exp(-im*(X*(Rfrom[2]-Rto[2]) - Y*(Rfrom[1]-Rto[1]))/q)
                end
            end
            # check minus in front of t here
            write(io, "      $(Rfrom[1])   $(Rfrom[2])   $n      $(Rto[1])   $(Rto[2])   $m        $(abs(t))      $(mod.(angle.(-1.0 .* t)/(2π) .+0.5 ,1.0))\n")
        end
    end
end

