include("wfunctions.jl")
using .Wannierize
include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil             

using Plots
using Measures
using Base.Threads
using ProgressMeter
using LinearAlgebra

outdir = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/wannier_out_plots/"


# physical inputs
p = 1
q = 2
a_nm = 5.0f0
U0 = 0.07f0

# numerical inputs
LLmax = 40
NXY = 48
Nuc = 30

# multiple usage values
phi = Float32(p/q)
a = a_nm * 1f-9
pi32 = Float32(π)

# grids

dx = q*a / (Nuc*q)
dy = a / Nuc

xrange_muc = Float32.(0:dx:(q*a - dx))
yrange_muc = Float32.(0:dy:(a - dy))
xrange = Float32.(0:dx:(2*q*a - dx))
yrange = Float32.(0:dy:(2*a - dy))


# test for 1/2; Q: how to reliably construct trial wavefunctions for abitrary k and q!

# trial wavefunction is k=0 Bloch in 1 unit cell
H0 = Hamil.get_full_ham(phi, 0f0*pi32, 0f0*pi32, U0, a, p, LLmax)
evals0, evecs0 = eigen(H0)
trialwf = Array{ComplexF32}(undef, Nuc*q, Nuc, q)
trialwf[:,:,1] = Wannierize.get_trial_from_u(1, xrange_muc,yrange_muc,0f0,0f0,evecs0[:,1],p,q,a,LLmax,dx,dy)
trialwf[:,:,2] = Wannierize.get_trial_from_u(2, xrange_muc,yrange_muc,0f0,0f0,evecs0[:,1],p,q,a,LLmax,dx,dy) 

# gaussian trial wf
# trialwf = Array{ComplexF32}(undef, Nuc*q, Nuc, q)
# for m in 1:q
#     normg = 0f0
#     for (i,x) in enumerate(xrange_muc)
#         for (j,y) in enumerate(yrange_muc)
#             g_xy = Wannierize.gaussian_LG(Float32(x),Float32(y),Float32(m-0.5)*a, 0f0,phi,U0,a)
#             g_xy *= exp((m-1)*im*y/a)
#             trialwf[i,j,m] = g_xy
#             normg += abs2(g_xy)
#         end
#     end
#     normg *= dx*dy
#     normg = sqrt(normg)
#     trialwf[:,:,m] ./= normg
# end

Xrange = range(0f0,2π*q,NXY*q+1)[1:end-1]
Yrange = range(0f0,2π,NXY+1)[1:end-1]

Rfrom = (0,0)
Rlist = [(Rx,Ry) for Rx=0:2 for Ry=0:2]

hopping_mats = zeros(ComplexF32, q, q, length(Rlist))


@showprogress for (iX,X) in enumerate(Xrange)
    for (jY,Y) in enumerate(Yrange)
        ky = X/(a*q) 
        kx = Y/(a*q)

        Umat, ens = Wannierize.get_U_matrix_and_energies(p,q,X,Y,U0,a,LLmax,Nuc,xrange_muc,yrange_muc,dx,dy,trialwf)

        Emat = diagm(ens[1:q])

        blochmat = Umat' * Emat * Umat

        for (i,Rto) in enumerate(Rlist)
            Rtox, Rtoy = Rto
            hopping_mats[:,:,i] .+= blochmat[:,:] .* exp(-im*a*(q*kx*(Rtox-Rfrom[1]) + ky*(Rtoy-Rfrom[2])))
        end

    end
end

hopping_mats[:,:,:] ./= (q*NXY^2)

outfile = joinpath(outdir, "t_mats_p$p-q$q-U0$U0-NXY$NXY-Nuc$Nuc-LLmax$LLmax.txt")

open(outfile, "w") do io
    for (i,Rto) in enumerate(Rlist)
        println(io, "Hopping $Rfrom --> $Rto")
        Wannierize.write_complex_matrix_in_io(io, hopping_mats[:,:,i])
        println(io,"\n")
    end
end








