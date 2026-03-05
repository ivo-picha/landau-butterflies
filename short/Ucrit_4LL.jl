# look at the critical potential for a Hofstadter butterfly gap to open in the 4-LL model, as a function of flux phi = 1/q
using Plots
using LinearAlgebra
using ProgressMeter
using Base.Threads
using Printf

include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil

#outdir = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/test"
outdir = "/users/ivoga/lh/out/test" # cluster path
mkpath(outdir)

U0list =  round.(collect(range(0.0010, 0.016, step = 0.0001)); digits=4)
a = 5f-9

qmax = 50
q_list = collect(range(2, qmax))

NXY = 64

energies = Array{Float32}(undef,length(U0list),length(q_list),4*NXY^2)

@showprogress @threads for U0 in U0list

    for q in q_list
        phi = 1f0/Float32(q)

        Xlist = collect(range(0f0, Float32(2π*q), length=(NXY+1)))[1:end-1]
        Ylist = collect(range(0f0, Float32(2π), length=(NXY+1)))[1:end-1]

        NLLmax = q+1
        enq = Float32[]
        for Y in Ylist
            for X in Xlist
                H = Hamil.get_full_ham(phi, X, Y, Float32(U0), a, 1, NLLmax)
                H = H[q-1:q+2, q-1:q+2] # get the 4x4 block around the filling 1 gap
                E = eigvals(Hermitian(H))
                append!(enq, E)
            end
        end
        energies[findfirst(==(U0), U0list), findfirst(==(q), q_list), :] = enq
    end
end

for (i,U0) in enumerate(U0list)
    plt = plot(framestyle=:box, xlabel="q", ylabel="E", size=(1200,800), title="U0 = $U0, a = 5 nm")
    for (j,q) in enumerate(q_list)
        energies_q = energies[i, j, :]
        scatter!(plt, q*ones(length(energies_q)), energies_q, markersize=2, color=:black, label="", markerstrokewidth=0)
    end
    savefig(plt, joinpath(outdir, "4LL-U0_$(@sprintf("%.6f", U0)).png"))
end