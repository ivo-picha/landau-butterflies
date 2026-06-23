# visualize the Wannier functions at different fluxes

# to do: make it work for any q 

using LinearAlgebra
using Plots
using ProgressMeter
using Base.Threads
using Measures

include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil                        
include(joinpath(dirname(@__DIR__),"funcs/states.jl"))
using .States 


p = 1
q = 2
a_nm = 5.0 # lattice constant in nm
U0 = 0.05f0 # potential strength

NXY = 251 # number of k-points in each direction per BZ UC
Ngrid = 60 # number of real space points in each lattice UC 
LLmax = 40 # max LL used


#output folder
output_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/wannier_out_plots/ph$p-$q-U$(abs(U0))-a$a_nm-LL$LLmax-NXY$NXY-Ngrid$Ngrid"
mkpath(output_folder)

a = Float32(a_nm*1f-9)
phi = Float32(p/q)

# ================================== comp

X_list = collect(range(0f0, Float32(2π*q), length = NXY*q+1))[1:end-1]
Y_list = collect(range(0f0, Float32(2π), length = NXY+1))[1:end-1]
XY_list = reshape(collect(Base.product(X_list, Y_list)),:)


xplotrange = range(-a, a*q+a, Ngrid*q) # maybe take a larger UC for better integral over "all space"
yplotrange = range(-a, 2*a, Ngrid)
dx = step(xplotrange)
dy = step(yplotrange)
xplotrange = Float32.(collect(xplotrange))
yplotrange = Float32.(collect(yplotrange))
# xyplotlist = reshape(collect(Base.product(xplotrange,yplotrange)),:)


# trial wavefunctions
println("Generating trial wavefunctions...")
g_array = Array{ComplexF32}(undef, Ngrid*q, Ngrid, q)
for m in 1:q
    normg = 0f0
    for (i,x) in enumerate(xplotrange)
        for (j,y) in enumerate(yplotrange)
            g_xy = States.gaussian_LG(x,y,Float32(m-1)*a + a/2, a/2,phi,U0,a)
            g_array[i,j,m] = g_xy
            normg += abs2(g_xy)
        end
    end
    normg *= dx*dy
    normg = sqrt(normg)
    g_array[:,:,m] ./= normg
end



pi32 = Float32(π)


wannier_array_re = [Threads.Atomic{Float32}(0.0f0) for i in 1:Ngrid*q, j in 1:Ngrid, k in 1:q]
wannier_array_im = [Threads.Atomic{Float32}(0.0f0) for i in 1:Ngrid*q, j in 1:Ngrid, k in 1:q]


println("Obtaining Wannier functions on $(Threads.nthreads()) thread(s)...")
@showprogress @threads for i = 1:NXY*q
    X = X_list[i]
    intermediate_wannier_array = zeros(ComplexF32, Ngrid*q, Ngrid, q)

    for j = 1:NXY
        Y = Y_list[j]

        states = Tuple{Float32, Vector{ComplexF32}}[] # (energy, eigenvec, X, Y)

        H = Hamil.get_full_ham(phi, X, Y, U0, a, p, LLmax)
        evals, evecs = eigen(H)
        for i in 1:q
            push!(states, (evals[i], evecs[:,i]))
        end

        # sort states by energy
        sort!(states, by = first)

        # momenta for convenience
        ky = X/(a*q)
        kx = Y/(a*q)

        # array to store wavefunctions
        wf_array = Array{ComplexF32}(undef, Ngrid*q, Ngrid, q)
        for m in 1:q
            norm = 0f0
            for (i,x) in enumerate(xplotrange)
                for (j,y) in enumerate(yplotrange)
                    wf_xy = States.get_eigenstate_XY(x,y,states[m],X,Y,p,q,a,LLmax) # full Bloch wf here
                    # stripping factor to reduce to periodic function in magnetic UC
                    peel = 1
                    #peel *= exp(-im*ky*y) #will be multiplied later anyways, 3 lines later
                    peel *= exp(im*kx*a*q*floor(x/(a*q))) #not needed?
                    peel *= exp(-im*2*pi32*(y/a)*p*floor(x/(a*q)))
                    wf_xy *= peel
                    # reconstruct new bloch function
                    wf_xy *= exp(im*kx*x)
                    wf_array[i,j,m] = wf_xy
                    norm += abs2(wf_xy)
                end
            end
            norm *= dx*dy
            norm = sqrt(norm)
            wf_array[:,:,m] ./= norm
        end


        # matrix Aₘₙ=⟨ψₘ|gₙ⟩
        Amat = Array{ComplexF32}(undef, q, q)
        for m in 1:q
            for n in 1:q
                Amn = sum(g_array[:,:,n] .* conj.(wf_array[:,:,m])) * dx*dy
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

        # if minimum(λ) < tol
        #     @warn "Overlap matrix nearly singular" minimum(λ) maximum(λ)
        # end

        inv_sqrt_λ = [x > tol ? 1 / sqrt(x) : 0.0 for x in λ]
        Smat_inv_sqrt = F.vectors * Diagonal(inv_sqrt_λ) * F.vectors'
        Umat = Amat * Smat_inv_sqrt

        # test
        Ucheck = norm(Umat'*Umat - I)
        if Ucheck > 1e-4
            println("Warning: U matrix is not unitary at X=$X, Y=$Y, norm: $Ucheck")
            # println(Umat'*Umat)
            # println("A matrix:")
            # println(round.(Amat;digits=6))
            # println("Energies:")
            # println(round.(energies;digits=8))
        end

        # add to FT
        for (i,x) in enumerate(xplotrange)
            for (j,y) in enumerate(yplotrange)
                for n = 1:q
                    Rx = Float32(n-1)*a + a/2
                    Ry = 0f0
                    w_xy_n = 1/(q*NXY^2)
                    w_xy_n *= exp(-im*(kx*Rx - ky*Ry))
                    w_xy_n *= sum([Umat[m,n] * wf_array[i,j,m] for m=1:q])

                    intermediate_wannier_array[i,j,n] += w_xy_n
                end
            end
        end
        
    end

    for (i,x) in enumerate(xplotrange)
        for (j,y) in enumerate(yplotrange)
            for n = 1:q
                iwa_xy_n = intermediate_wannier_array[i,j,n]
                Threads.atomic_add!(wannier_array_re[i,j,n], real(iwa_xy_n))
                Threads.atomic_add!(wannier_array_im[i,j,n], imag(iwa_xy_n))
            end
        end
    end
end

g_array = nothing;


# normalize wannier

wannier_rebuilt = Array{ComplexF32}(undef, Ngrid*q, Ngrid, q)
for (i,x) in enumerate(xplotrange)
    for (j,y) in enumerate(yplotrange)
        for n = 1:q
            wannier_rebuilt[i,j,n] = wannier_array_re[i,j,n][] + im*wannier_array_im[i,j,n][]
        end
    end
end

plots_vec = Plots.Plot[];

for n = 1:q
    absdat = abs.(wannier_rebuilt[:,:,n])' / maximum(abs.(wannier_rebuilt[:,:,n]))
    plt_abs = heatmap(xplotrange./a, yplotrange./a, absdat, 
        xaxis = "x/a", yaxis = "y/a", title = "|w|", color = :viridis)
    plt_arg = heatmap(xplotrange./a, yplotrange./a, angle.(wannier_rebuilt[:,:,n])', 
        xaxis = "x/a", yaxis = "y/a", title = "arg(w)", color = :hsv)

    push!(plots_vec, plt_abs)
    push!(plots_vec, plt_arg)
end


plt_combi = plot(plots_vec..., layout = (q,2), dpi = 120, size = (900,450), measures = 5mm)

savefig(plt_combi, joinpath(output_folder, "test2.png"))

