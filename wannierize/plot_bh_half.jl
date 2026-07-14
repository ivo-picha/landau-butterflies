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
LLmax = 40
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

Xrange = range(0f0,2π*q,NXY*q+1)[1:end-1]
Yrange = range(0f0,2π,NXY+1)[1:end-1]


bloch11_array = Array{ComplexF32}(undef, NXY*q, NXY)
bloch12_array = Array{ComplexF32}(undef, NXY*q, NXY)
bloch21_array = Array{ComplexF32}(undef, NXY*q, NXY)
bloch22_array = Array{ComplexF32}(undef, NXY*q, NXY)



@showprogress for (iX,X) in enumerate(Xrange)
    for (jY,Y) in enumerate(Yrange)
        ky = X/(a*q) 
        kx = Y/(a*q)

        Umat, ens = Wannierize.get_U_matrix_and_energies(p,q,X,Y,U0,a,LLmax,Nuc,xrange_muc,yrange_muc,dx,dy,trialwf)

        Emat = diagm(ens[1:q])

        blochmat = Umat' * Emat * Umat

        bloch11_array[iX,jY] = blochmat[1,1]
        bloch12_array[iX,jY] = blochmat[1,2]
        bloch21_array[iX,jY] = blochmat[2,1]
        bloch22_array[iX,jY] = blochmat[2,2]

    end
end

plt_b11_abs = heatmap(Xrange, Yrange, abs.(bloch11_array[:,:])', color=:viridis, aspect_ratio=1, xaxis = "X", yaxis = "Y")
plt_b11_arg = heatmap(Xrange, Yrange, angle.(bloch11_array[:,:])', color=:hsv, aspect_ratio=1, xaxis = "X", yaxis = "Y")
pltcombi_b11 = plot(plt_b11_abs, plt_b11_arg, layout = (2,1), size = (q*500,1100), framestyle=:box)

plt_b12_abs = heatmap(Xrange, Yrange, abs.(bloch12_array[:,:])', color=:viridis, aspect_ratio=1, xaxis = "X", yaxis = "Y")
plt_b12_arg = heatmap(Xrange, Yrange, angle.(bloch12_array[:,:])', color=:hsv, aspect_ratio=1, xaxis = "X", yaxis = "Y")
pltcombi_b12 = plot(plt_b12_abs, plt_b12_arg, layout = (2,1), size = (q*500,1100), framestyle=:box)

plt_b21_abs = heatmap(Xrange, Yrange, abs.(bloch21_array[:,:])', color=:viridis, aspect_ratio=1, xaxis = "X", yaxis = "Y")
plt_b21_arg = heatmap(Xrange, Yrange, angle.(bloch21_array[:,:])', color=:hsv, aspect_ratio=1, xaxis = "X", yaxis = "Y")
pltcombi_b21 = plot(plt_b21_abs, plt_b21_arg, layout = (2,1), size = (q*500,1100), framestyle=:box)

plt_b22_abs = heatmap(Xrange, Yrange, abs.(bloch22_array[:,:])', color=:viridis, aspect_ratio=1, xaxis = "X", yaxis = "Y")
plt_b22_arg = heatmap(Xrange, Yrange, angle.(bloch22_array[:,:])', color=:hsv, aspect_ratio=1, xaxis = "X", yaxis = "Y")
pltcombi_b22 = plot(plt_b22_abs, plt_b22_arg, layout = (2,1), size = (q*500,1100), framestyle=:box)

plt_all_bloch = plot(pltcombi_b11, pltcombi_b12, pltcombi_b21, pltcombi_b22, layout = (2,2), size = (q*500,1200), margins = 3mm)


#plt_comb = plot(plt_wannier1,plt_wannier2,layout=(1,2), size=(1000,500))

#savefig(plt_comb, joinpath(outdir,"test_wannier_plusY_p$p-q$q-U0$U0-LLmax$LLmax-NXY$NXY-Nuc$Nuc.png"))






