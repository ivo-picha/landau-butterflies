using LinearAlgebra
using Plots
using ProgressMeter
using KernelDensity

# params in meV; U0 is 50meV
mu = -30.0846
tNN = 0.615
tNNN = 0.018
tdNNN = 0.026

# k-space resolution
Nkxy = 128

# kx and ky dimensionless from 0 to 2pi
function get_h_tb(kx, ky, mu, tNN, tNNN, tdNNN)
    h = zeros(ComplexF64, 2, 2)
    h[1,1] = mu/2 + exp(im*ky)*tNN*exp(im*π/2) - exp(im*kx)*tNNN - exp(im*ky*2)*tNNN
    h[2,2] = mu/2 + exp(im*ky)*tNN*exp(-im*π/2) - exp(im*kx)*tNNN - exp(im*ky*2)*tNNN
    h[1,2] = - exp(im*kx)*tNN - exp(im*kx)*exp(im*ky)*tdNNN - exp(im*ky)*tdNNN         # controversial minus in front
    h[2,1] = tNN - exp(-im*kx)*exp(im*ky)*tdNNN - exp(im*ky)*tdNNN
    return Hermitian(h+h')
end

function get_DOS_x_d(data; res=100)
    eta = maximum(data)-minimum(data)
    eta /= res
    kd = kde(data, bandwidth = eta)
    return kd.x, kd.density
end

kxy_range = range(0,2π,Nkxy+1)[1:end-1]

energies = Float64[];

@showprogress for kx in kxy_range
    for ky in kxy_range
        append!(energies, eigvals(get_h_tb(kx, ky, mu, tNN, tNNN, tdNNN)))
    end
end

sort!(energies)
en_x, en_d = get_DOS_x_d(energies)

# plt = Plots.plot(
#     framestyle = :box,
#     xaxis = "Energy [meV]",
#     yaxis = "DOS",
#     title = "ϕ = 1/2, U₀ = 50 meV, a = 5 nm"

# )

plot!(plt, en_x,en_d, label = "NN+dNNN+NNN", lw = 2, color = 3, legend = :top)



## add original spectrum
include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil

energies_original = Float32[];
@showprogress for X in kxy_range
    X *= 2
    for Y in kxy_range
        H = Hamil.get_full_ham(phi, Float32(X), Float32(Y), 0.05f0, 5.0f0, 1, 40)
        append!(energies_original, eigvals(H)[1:2])
    end
end

sort!(energies_original)
en_o_x, en_o_d = get_DOS_x_d(energies_original)

plot!(plt, en_x,en_d, label = "original", color = :black, style = :dot, lw = 4)

