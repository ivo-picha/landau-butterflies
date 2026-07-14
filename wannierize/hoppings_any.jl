include("wfunctions.jl")
using .Wannierize
include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil             

# using Plots
# using Measures
using Base.Threads
using ProgressMeter
using LinearAlgebra

outdir = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/wannier_out_hops/"
#outdir = "/users/ivoga/lh/out/wannier_out" # cluster path

mkpath(outdir)

args = ARGS
# p q a_nm U0 LLmax NXY Nuc
if length(args) != 7
    println("ERROR: Arguments should be in order < p q a_nm U0 LLmax NXY Nuc >.")
end

p = parse(Int, args[1])
q = parse(Int, args[2])
a_nm = parse(Float32, args[3])
U0 = parse(Float32, args[4])
LLmax = parse(Int, args[5])
NXY = parse(Int, args[6])
Nuc = parse(Int, args[7])


# # physical inputs
# p = 1
# q = 3
# a_nm = 5.0f0
# U0 = 0.05f0

# # numerical inputs
# LLmax = 60
# NXY = 64
# Nuc = 30

# multiple usage values
phi = Float32(p/q)
a = a_nm * 1f-9
pi32 = Float32(π)

# grids
xrange_muc = range(0f0,q*a,q*Nuc)
yrange_muc = range(0f0,a,Nuc)
xrange = range(0f0,2*q*a,q*Nuc)
yrange = range(0f0,2*a,Nuc)
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


@showprogress @threads for X in Xrange
    for Y in Yrange
        ky = X/(a*q) 
        kx = Y/(a*q)

        Umat, ens = Wannierize.get_U_matrix_and_energies(p,q,X,Y,U0,a,LLmax,Nuc,xrange_muc,yrange_muc,dx,dy,trial_wfs)

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
