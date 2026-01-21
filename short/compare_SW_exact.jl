# compare exact Hofstadter bands with Schrieffer-Wolff approximation

include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis
using Plots
using LinearAlgebra


r = 0.001 # ratio U0/ΔE

a = 5.0e-9
p = 2
q = 1
NLL_max = 10
phi = p/q
dE = Hamil.E_LL(1,Float32(phi),Float32(a)) - Hamil.E_LL(0,Float32(phi),Float32(a))
U0 = r*dE

# set up path in BZ
Nk = 200;
kpath = Tuple{Float64,Float64}[];
# add center -> edge1
for j=0:(Nk-1)
    kj = (0.0,π *j/Nk)
    push!(kpath, kj)
end
# add edge1 -> corner
for j=0:(Nk-1)
    kj = (0.0,π) .+ ((π,π).-(0.0,π)).*(j/Nk) 
    push!(kpath, kj)
end
# add corner -> gamma
for j=0:(Nk-1)
    kj = (π,π) .+ ((0.0,0.0).-(π,π)).*(j/Nk) 
    push!(kpath, kj)
end
# add gamma -> edge2
for j=0:(Nk-1)
    kj = (0.0,0.0) .+ ((π,0.0).-(0.0,0.0)).*(j/Nk) 
    push!(kpath, kj)
end



# get list of energies
ens = [];
for (k1,k2) in kpath
    H = Hamil.get_full_ham(Float32(phi), Float32(k1), Float32(k2), Float32(U0), Float32(a), p, NLL_max)
    push!(ens, sort(eigvals(H)))
end
bandes = hcat(ens...)


plt = Plots.plot(tickfontsize = 15,   # increase tick label size
    guidefontsize = 14,  # axis label size
    legendfontsize = 11, # legend text size
    titlefontsize = 15,   # title size
    size = (350,400)
)
xticks!(plt,[0,Nk,2*Nk,3*Nk,4*Nk],["Γ","Xʸ","M","Γ","Xˣ"])
yticks!(plt,[0],[""])
ylabel!(plt, "Energy")
for j = 1:p
    plot!(plt, bandes[j,:], label = "B$j, exact", framestyle = :box, color = j,lw = 2, legend=:right)
end




# add approximate Schrieffer Wolf solution, up to first order in U0/ΔE
function HSW(k1,k2)
    th00 = Float64(Hamil.T(0,0,Float32(phi)))
    th10 = Float64(Hamil.T(1,0,Float32(phi)))
    En0 = Hamil.E_LL(0,Float32(phi),Float32(a))
    En1 = Hamil.E_LL(1,Float32(phi),Float32(a))
    dE = En1 - En0
    h11 = En0+U0*th00*cos(k1/2)-U0^2*th10^2/dE*(sin(k1/2)^2+sin(k2/2)^2)
    h22 = En0-U0*th00*cos(k1/2)-U0^2*th10^2/dE*(sin(k1/2)^2+sin(k2/2)^2)
    h12 = U0*th00*cos(k2/2) + im* 2*U0^2*th10^2/dE*sin(k1/2)*sin(k2/2)
    return Hermitian([h11 h12; conj(h12) h22])
end


# get list of energies
ensSW = [];
for (k1,k2) in kpath
    H = HSW(k1,k2)
    push!(ensSW, sort(eigvals(H)))
end
bandesSW = hcat(ensSW...)

for j = 1:p
    plot!(plt, bandesSW[j,:], label = "B$j, S-W", framestyle = :box, linestyle = :dot, color = j, lw = 2)
end
title!(plt, "U₀/Δε = $r")

#plt
savefig(plt, joinpath(dirname(@__DIR__),"/home/ivoga/Documents/PhD/Landau_Hofstadter/figs/Schrieffer-Wolff-phi2/bands_r$r.png"))

