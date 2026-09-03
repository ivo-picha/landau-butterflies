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
#outdir = "/users/ivoga/lh/out/wannier_out" # cluster path

mkpath(outdir)

# args = ARGS
# p q a_nm U0 LLmax NXY Nuc
# if length(args) != 7
#     println("ERROR: Arguments should be in order < p q a_nm U0 LLmax NXY Nuc >.")
# end

# p = parse(Int, args[1])
# q = parse(Int, args[2])
# a_nm = parse(Float32, args[3])
# U0 = parse(Float32, args[4])
# LLmax = parse(Int, args[5])
# NXY = parse(Int, args[6])
# Nuc = parse(Int, args[7])


# physical inputs
p = 3
q = 2
a_nm = 5.0f0
U0 = 0.05f0

# numerical inputs
LLmax = 40
NXY = 16
Nuc = 60

# multiple usage values
phi = Float32(p/q)
a = a_nm * 1f-9
pi32 = Float32(π)

# grids
xrange_muc = range(0f0,q*a,length=q*Nuc)
yrange_muc = range(0f0,a,length=Nuc)
xrange = range(0f0,2*q*a,length=2*q*Nuc)
yrange = range(0f0,2*a,length=2*Nuc)
dx = step(xrange_muc)
dy = step(yrange_muc)


# trial wfunctions

trial_wfs = Wannierize.get_all_trial_wavefunctions(p,q,U0,a,LLmax,xrange_muc,yrange_muc,dx,dy)
#Wannierize.plot_array(xrange_muc,yrange_muc, trial_wfs[:,:,5],a,q)


Xrange = range(0f0,2π*q,NXY*q+1)[1:end-1]
Yrange = range(0f0,2π,NXY+1)[1:end-1]

Rfrom = (0,0)
Rlist = [(Rx,Ry) for Rx=0:2 for Ry=0:2]

hopping_mats = zeros(ComplexF32, q, q, length(Rlist))
wannier_wf_array = zeros(ComplexF32, length(xrange), length(yrange), q)

@showprogress @threads for X in Xrange
    for Y in Yrange
        ky = X/(a*q) 
        kx = Y/(a*q)

        la = Wannierize.get_loewdin_array(p,q,X,Y,U0,a,LLmax,Nuc,xrange,yrange,xrange_muc,yrange_muc,dx,dy,trial_wfs)

        for m = 1:q
            ft_factor = exp(-im*(kx*(m-0.5)*a + ky*(0)*a)) / (NXY^2 *q)
            wannier_wf_array[:,:,m] .+= ft_factor.*la[:,:,m]
        end

    end
end

##
plt2 = Wannierize.plot_array(xrange_muc,yrange_muc,wannier_wf_array[1:Nuc*q,1:Nuc,1],a,q)
yaxis!(plt2, "y/a")



# plt_list = [];
# for m = 1:q
#     push!(plt_list, Wannierize.plot_array(xrange,yrange,wannier_wf_array[:,:,m],a,q))
# end

plt_comb = plot(plt1, layout=(1,1), size=(320,600), dpi = 200)

savefig(plt_comb, joinpath(outdir,"plot_wannier_p$p-q$q-U0$U0-LLmax$LLmax-NXY$NXY-Nuc$Nuc.pdf"))
