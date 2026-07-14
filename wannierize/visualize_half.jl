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
U0 = 0.05f0

# numerical inputs
LLmax = 8
NXY = 32
Nuc = 20

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
H0 = Hamil.get_full_ham(phi, 0f0, 0f0, U0, a, p, LLmax)
evals0, evecs0 = eigen(H0)
trialwf = Array{ComplexF32}(undef, Nuc*q, Nuc, q)
trialwf[:,:,1] = Wannierize.get_trial_from_u(1, xrange_muc,yrange_muc,0f0,0f0,evecs0[:,1],p,q,a,LLmax,dx,dy)
trialwf[:,:,2] = Wannierize.get_trial_from_u(2, xrange_muc,yrange_muc,0f0,0f0,evecs0[:,1],p,q,a,LLmax,dx,dy)

#Wannierize.plot_array(xrange_muc, yrange_muc, trialwf[:,:,2],a,q)

# gaussian trial wf
# trialwf = Array{ComplexF32}(undef, Nuc*q, Nuc, q)
# for m in 1:q
#     normg = 0f0
#     for (i,x) in enumerate(xrange_muc)
#         for (j,y) in enumerate(yrange_muc)
#             g_xy = Wannierize.gaussian_LG(x,y,Float32(m-0.5)*a, 0f0,phi,U0,a)
#             trialwf[i,j,m] = g_xy
#             normg += abs2(g_xy)
#         end
#     end
#     normg *= dx*dy
#     normg = sqrt(normg)
#     trialwf[:,:,m] ./= normg
# end

Xrange = range(0f0,2π*q,NXY+1)[1:end-1]
Yrange = range(0f0,2π,NXY+1)[1:end-1]

wannier_wf_array = zeros(ComplexF32, 2*Nuc*q, 2*Nuc, q)

@showprogress for X in Xrange
    for Y in Yrange
        ky = X/(a*q) 
        kx = Y/(a*q)

        la = Wannierize.get_loewdin_array(p,q,X,Y,U0,a,LLmax,Nuc,xrange,yrange,xrange_muc,yrange_muc,dx,dy,trialwf)

        for m = 1:q
            ft_factor = exp(-im*(kx*(m-0.5)*a + ky*(0)*a)) / (NXY^2 *q)
            wannier_wf_array[:,:,m] .+= ft_factor.*la[:,:,m]
        end

    end
end


plt_wannier1 = Wannierize.plot_array(xrange,yrange,wannier_wf_array[:,:,1],a,q)
plt_wannier2 = Wannierize.plot_array(xrange,yrange,wannier_wf_array[:,:,2],a,q)
plt_comb = plot(plt_wannier1,plt_wannier2,layout=(1,2), size=(1000,500))

savefig(plt_comb, joinpath(outdir,"test_wannier_p$p-q$q-U0$U0-LLmax$LLmax-NXY$NXY-Nuc$Nuc.png"))






