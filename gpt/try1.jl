using Plots
using ProgressMeter
using LinearAlgebra

include("FluxTightBinding.jl")
using .FluxTightBinding


Nk = 64
qmax = 64
p_range = 1:1:2*qmax

energies = Float64[];
phis = Float64[];
@showprogress for pn in p_range
    gcdpq = gcd(pn,qmax)
    p = div(pn,gcdpq)
    q = div(qmax,gcdpq)
    phi = Float64(p/q)

    kx_range = range(-π/q,π/q, Nk)
    ky_range = range(-π,π, Nk)

    energies_at_phi = Float64[];
    for kx in kx_range
        for ky in ky_range
            h = FluxTightBinding.tb_hamiltonian(p,q,kx,ky)
            append!(energies_at_phi, eigvals(h))
        end
    end
    append!(energies, energies_at_phi)
    append!(phis,[phi for j=1:length(energies_at_phi)])
end


# plot
plt = Plots.scatter(phis,energies,
        xlabel = "ϕ = p/q", ylabel = "Energy [meV]", title = "",
        label = "", ms = 0.5, size = (800,600), color = :red, markerstrokewidth=0,
        xaxis = false, yaxis = false, grid = false, background_color = :transparent
);
out_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/modified_butterflies/"
mkpath(out_folder)
savefig(plt,joinpath(out_folder,"modb_v3_from0to2_U0-50meV.png"))


## params

alphas = range(0.0,3,100)
tdNNNlist = Float64[]
for alpha in alphas 
   append!(tdNNNlist, FluxTightBinding.tb_parameters(alpha).tdNNN) 
end
plot(alphas,tdNNNlist)


## test

scatter([
    (2//7,  0.927),
    (3//7,  0.865),
    (1//2,  0.001),  # q=2: orientation-aliased
    (3//5,  0.766),
    (5//8,  0.817),
    (5//7,  missing), # failed
    (3//4,  0.562),
    (6//7,  0.799),
    (1//1, -0.987),
    (7//6, -0.800),
    (7//5,  missing), # failed
    (5//3, -0.130),
    (7//4, -0.064),
    (2//1, -0.010),
    (8//3,  0.942),
    (3//1,  0.169),   # near zero: phase unstable
], xaxis = "ϕ", yaxis = "arg(t_diag)/π", label = "")

scatter([
    (2//7, 57.5),
    (3//7, 58.2),
    (1//2, 25.7),  # q=2: orientation-aliased
    (3//5, 70.1),
    (5//8, 39.2),
    (5//7, missing), # failed
    (3//4, 72.4),
    (6//7, 18.6),
    (1//1, 25.0),
    (7//6, 14.1),
    (7//5, missing), # failed
    (5//3, 18.1),
    (7//4, 17.5),
    (2//1, 13.0),
    (8//3,  5.0),
    (3//1,  2.0),   # near origin
], xaxis = "ϕ", yaxis = "|t_diag|", label = "")

## my attempt

function print_diagonal_phases_pi(p,q)
    phi(m) = round((mod(2π * (p/q) * (m+1) + π,2π) - π)/π ; digits = 6)
    println("Diagonal phases for ϕ = $p/$q:")
    for m=0:(q-1) 
        println("m = $m: $(phi(m))")
    end
end

print_diagonal_phases_pi(3,2)

dNNN_results = [
    (0, 1, [0],                   false, 1.474e-6,  0.000, 1.00),
    (1, 7, [6, 1, 5, 0, 4, 2, 3], true,  4.125e-6,  0.037, 0.84),
    (1, 6, [0, 5, 1, 4, 2, 3],    true,  4.875e-6, -0.025, 0.99),
    (1, 5, [0, 4, 3, 1, 2],       true,  7.327e-6, -0.020, 0.88),
    (1, 4, [0, 3, 1, 2],          false, 8.882e-6, -0.027, 0.99),
    (2, 7, [3, 6, 0, 2, 4, 1, 5], false, 1.032e-4,  0.008, 0.52),
    (1, 3, [0, 2, 1],             false, 1.432e-5, -0.029, 0.99),
    (2, 5, [2, 0, 4, 1, 3],       true,  1.004e-4, -0.286, 0.38),
    (3, 7, [2, 4, 0, 6, 5, 1, 3], true,  6.013e-5, -0.246, 0.63),
    (1, 2, [1, 0],                false, 2.584e-5,  0.031, 1.00),
    (4, 7, [3, 5, 1, 6, 0, 4, 2], false, 6.500e-5, -0.082, 0.91),
    (3, 5, [1, 3, 4, 0, 2],       false, 8.849e-5, -0.088, 0.80),
    (2, 3, [1, 0, 2],             true,  7.798e-5, -0.150, 0.58),
    (5, 7, [5, 1, 4, 2, 0, 3, 6], true,  1.834e-5, -0.279, 0.73),
    (3, 4, [1, 2, 0, 3],          false, 9.727e-5, -0.207, 0.49),
    (5, 6, [2, 3, 0, 4, 5, 1],    true,  4.262e-5,  0.083, 0.51),
    (6, 7, [3, 2, 4, 0, 1, 5, 6], true,  4.356e-5, -0.054, 0.61),
    (1, 1, [0],                   false, 2.500e-5,  0.000, 1.00),
    (7, 6, [3, 2, 4, 0, 1, 5],    true,  1.507e-5,  0.251, 0.85),
    (4, 3, [1, 2, 0],             false, 2.563e-5,  0.512, 1.00),
    (7, 5, [1, 3, 0, 2, 4],       true,  3.814e-5,  0.846, 0.45),
    (3, 2, [1, 0],                false, 3.877e-5,  0.495+0.2, 1.00),
    (8, 5, [2, 4, 0, 1, 3],       false, 1.465e-5,  0.907, 0.99),
    (5, 3, [0, 2, 1],             false, 1.548e-5,  0.926, 0.96),
    (7, 4, [3, 0, 2, 1],          false, 1.557e-5,  0.967, 0.98),
    (2, 1, [0],                   false, 1.300e-5,  1.000, 1.00),
    (7, 3, [0, 2, 1],             false, 1.491e-6,  0.082, 0.13),
    (5, 2, [0, 1],                false, 1.949e-6, -0.150, 0.55),
    (8, 3, [1, 0, 2],             false, 3.659e-6, -0.076, 0.98),
    (3, 1, [0],                   false, 2.000e-6, -0.829, 1.00),
]

plt_my1abs=plot(framestyle=:box,xaxis="ϕ=p/q",yaxis="|t_diag|")
plt_my1arg=plot(framestyle=:box,xaxis="ϕ=p/q",yaxis="arg(t_diag)")
for dnnn in dNNN_results
    phi = dnnn[1]/dnnn[2]
    yarg = dnnn[6]
    yabs = dnnn[5]
    if dnnn[end]>0.7 #|| !dnnn[4]
        scatter!(plt_my1abs,(phi,yabs),label="",color=1,marker=:diamond, ms=6)
        scatter!(plt_my1arg,(phi,yarg),label="",color=1,marker=:diamond, ms=6)
    end
end
plt_my1comb = plot(plt_my1abs,plt_my1arg,layout=(2,1),size=(600,500),dpi=120)