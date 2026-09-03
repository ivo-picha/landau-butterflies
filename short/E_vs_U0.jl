# plot energy vs U0 for a fixed flux 

using Plots
using LinearAlgebra
using ProgressMeter
using Measures

include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis


p = 1
q = 2
a = 5f-9

# U0 in meV
U0_i = 0.0f0
U0_f = 22.0f0

LLmax = 50
LLshow = 6
NXY = 64
NU = 50

U0range = range(U0_i,U0_f,NU)
Xrange = range(0f0,Float32(2π*q),NXY)
Yrange = range(0f0,Float32(2π),NXY)
phi = Float32(p/q)

energies = Float32[];
U0s = Float32[];
@showprogress for U0 in U0range
    for X in Xrange
        for Y in Yrange
            H = Hamil.get_full_ham(phi, X, Y, U0/1000, a, p, LLmax)
            evH = eigvals(H)[1:LLshow]
            append!(energies,evH)
            append!(U0s,[U0 for i in eachindex(evH)])
        end
    end
end

## 2d plot
plt = Plots.scatter(U0s,energies.*1000,
    framestyle = :box,
    xaxis = "U₀ [meV]",
    title = "ϕ = $p/$q",
    yaxis = "Energy [meV]",
    markersize = 1,
    color = :black,
    markerstrokewidth = 0,
    size = (800,500),
    dpi = 300,
    label = "",
);

figpath = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/E_vs_U0_at_phi_$(round(phi;digits=2)).png"
savefig(plt, figpath)

## DOS 3d plot
using KernelDensity
using CairoMakie

fig = CairoMakie.Figure(size = (950, 700))

ax = CairoMakie.Axis3(
    fig[1, 1],
    title = "ϕ = $p/$q",
    titlesize = 30,
    xlabel = "U₀ [meV]",
    ylabel = "Energy [meV]",
    zlabel = "DOS",
    zticklabelsvisible = false,
    xlabelsize = 30,
    xticklabelsize = 25,
    ylabelsize = 30,
    yticklabelsize = 25,
    zlabelsize = 30,
    azimuth = 0.04*π,
    elevation = π / 4,
    aspect = (1.2, 2, 1)
)

for (j,U0) in enumerate(U0range)
    en_j = energies[(j-1)*p*LLshow*NXY^2 + 1 : j*p*LLshow*NXY^2]
    sort!(en_j)
    eta_j = (maximum(en_j)-minimum(en_j))/300
    kd_j = kde(en_j, bandwidth = eta_j)

    CairoMakie.lines!(
        ax,
        fill(U0,length(kd_j.x)),
        (kd_j.x).*1000,
        kd_j.density,
        color = :black,
        linewidth = 2
    )

end

fig
save("dos_E_U0_phi_1_2.png", fig, px_per_unit = 2)

## DOS 2d plot

plt_dos2d = Plots.plot(
    framestyle = :box,
    ylabel = "← U₀     DOS",
    xlabel = "Energy [meV]",
    yticks = false,
    tickfontsize = 12,
    margins = 5mm,
    xlims = (-7,40)
    )

for (j,U0) in enumerate(U0range)

    #iseven(j)&&continue

    en_j = energies[(j-1)*p*LLshow*NXY^2 + 1 : j*p*LLshow*NXY^2]
    sort!(en_j)
    eta_j = (maximum(en_j)-minimum(en_j))/200
    kd_j = kde(en_j, bandwidth = eta_j)

    Plots.plot!(
        plt_dos2d,
        (kd_j.x).*1000,
        kd_j.density .- (U0*70),
        color = :black,
        linewidth = 2,
        label = ""
    )

end

plt_dos2d