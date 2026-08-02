using Plots
using LinearAlgebra
using Interpolations
using ProgressMeter

#include("extended_hofstadter.jl")

# params
Nk = 8 # num k-points
qmax = 120

# fluxes
pqlist = [(0,1),
          (1, 8), (1, 7), (1, 6), (1, 5), (1, 4), (2, 7), (1, 3), (3, 8), (2, 5), (3, 7), 
          (1, 2), (4, 7), (3, 5), (5, 8), (2, 3), (5, 7), (3, 4), (4, 5), (5, 6), (6, 7),
          (7, 8), (1, 1), (8, 7), (7, 6), (6, 5), (5, 4), (4, 3), (7, 5), (3, 2), (8, 5),
          (5, 3), (7, 4), (2, 1), (7, 3), (5, 2), (8, 3), (3, 1)]
philist = [p//q for (p,q) in pqlist]

# data 
# chemical potential
mu_list = [-0.0701928,
       missing, missing, -0.07004, -0.0703, -0.07014, missing, -0.07011, missing, -0.070094, missing,
       -0.070028, missing, -0.069975, missing, -0.069925, missing, -0.06986, missing, missing, missing,
       missing, -0.069422, missing, missing, -0.06936, missing, -0.069156, missing, -0.068867, missing,
       missing, -0.068451, -0.06779, -0.067127, -0.066672, -0.066205, -0.065098]
# NN hopping
tnn_list = [0.0002035504457158282,
        missing, missing, missing, 0.000241, 0.000248, missing, 0.000246, missing, 0.000244, missing,
        0.000241, missing, 0.000236, missing, 0.000235, missing, 0.000234, missing, missing, missing,
        missing, 0.000201, missing, missing, missing, missing, 0.000194, missing, 0.000185, missing,
        missing, 0.000166, 0.000164, 0.000125, 0.000113, 0.0001, 7.9e-5]
# NNN hopping
tnnn_list = [4.32333950e-5,
         missing, missing, missing, 4.0e-6, 5.0e-6, missing, 2.0e-6, missing, 2.0e-6, missing,
         2.0e-6, missing, 8.0e-6, missing, 2.0e-6, missing, 5.0e-6, missing, missing, missing,
         missing, 5.2e-5, missing, missing, missing, missing, 4.0e-6, missing, 7.0e-6, missing,
         missing, 1.0e-6, 1.7e-5, 1.0e-6, 1.6e-5, 5.0e-6, 2.3e-5] 
# d-NNN hopping
dNNN_results_80meV = [
    (0, 1, [0],                    false, 2.8144e-07, +0.0000, 1.000),
    (1, 7, [5, 4, 2, 6, 1, 3, 0], true,  1.2004e-05, -0.0216, 0.455),
    (1, 6, [0, 4, 1, 3, 5, 2],    true,  1.9266e-05, +0.2445, 0.850),
    (1, 5, [1, 0, 4, 2, 3],       true,  2.3752e-06, -0.3035, 0.953),
    (1, 4, [0, 3, 1, 2],          false, 9.9431e-07, -0.0186, 0.995),
    (1, 3, [0, 2, 1],             false, 1.9990e-06, -0.0214, 0.995),
    (2, 5, [2, 4, 0, 1, 3],       false, 3.8455e-05, +0.0083, 0.489),
    (1, 2, [1, 0],                false, 2.9899e-06, +0.0208, 0.998),
    (4, 7, [3, 5, 1, 6, 0, 4, 2], false, 1.5043e-05, -0.0427, 0.860),
    (3, 5, [1, 3, 0, 4, 2],       false, 2.2686e-05, -0.0525, 0.715),
    (2, 3, [1, 0, 2],             false, 4.7757e-05, -0.1930, 0.620),
    (5, 7, [1, 5, 4, 2, 6, 0, 3], false, 1.4533e-05, -0.1314, 0.456),
    (3, 4, [2, 1, 0, 3],          false, 3.0432e-05, -0.3167, 0.320),
    (4, 5, [2, 4, 1, 3, 0],       true,  8.4261e-05, -0.4679, 0.875),
    (1, 1, [0],                    false, 1.4000e-05, +0.9637-1, 1.000),
    (8, 7, [3, 5, 4, 0, 6, 2, 1], true,  1.2384e-04, -0.8272+1, 0.843),
    (7, 6, [3, 2, 4, 0, 1, 5],    true,  2.9066e-06, +0.1092, 0.786),
    (6, 5, [2, 1, 3, 4, 0],       true,  1.1620e-05, +0.5469, 0.938),
    (5, 4, [2, 1, 3, 0],          false, 7.6069e-06, +0.4579, 0.963),
    (4, 3, [1, 2, 0],             false, 1.0254e-05, +0.4810, 0.992),
    (7, 5, [1, 3, 4, 0, 2],       false, 3.6680e-06, +0.5371, 0.897),
    (3, 2, [0, 1],                false, 1.4493e-05, +0.4739, 0.998),
    (8, 5, [2, 0, 4, 1, 3],       false, 2.5917e-06, +0.9004, 0.953),
    (5, 3, [0, 1, 2],             true,  1.8888e-05, -0.3162+1, 0.941),
    (7, 4, [0, 3, 2, 1],          false, 3.6600e-06, -0.9539+2, 0.872),
    (2, 1, [0],                    false, 6.0000e-06, -0.6863+1.6, 1.000),
    (7, 3, [0, 2, 1],             false, 9.9946e-07, +0.4440, 0.768),
    (5, 2, [1, 0],                false, 9.9979e-07, +0.5562, 0.978),
    (8, 3, [1, 0, 2],             false, 8.2848e-07, -0.1909, 1.000),
    (3, 1, [0],                    true,  2.7000e-05, -0.5586, 1.000),
]


dNNN_results_50meV = [
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
# Each tuple contains:
# (
#   p,
#   q,
#   index_to_m,
#   mixed,
#   mean_abs_t_dNNN,
#   residual_phase_over_pi,
#   coherence
# )
#
# residual_phase = arg(t_dNNN) - phase_Hofstadter - arg(t_xNN)
# Actual residual phase in radians is π * residual_phase_over_pi.

dNNN_results_20meV = [
    (0, 1, [0],                    false, 1.4311e-05, +0.0000, 1.000),
    (1, 7, [6, 5, 1, 0, 4, 2, 3], true,  5.3041e-05, +0.0227, 0.986),
    (1, 6, [1, 5, 4, 0, 2, 3],    true,  1.5747e-04, -0.1017, 0.739),
    (1, 5, [0, 4, 3, 1, 2],       true,  9.3446e-05, -0.0314, 0.983),
    (1, 4, [0, 3, 1, 2],          false, 1.3184e-04, -0.0349, 0.984),
    (2, 7, [3, 6, 0, 2, 4, 5, 1], false, 4.1450e-04, +0.0096, 0.725),
    (1, 3, [0, 2, 1],             false, 1.9721e-04, -0.0387, 0.987),
    (2, 5, [2, 4, 0, 1, 3],       true,  3.5530e-04, -0.1104, 0.622),
    (3, 7, [4, 2, 0, 6, 5, 3, 1], true,  4.7034e-04, -0.0774, 0.925),
    (1, 2, [1, 0],                false, 3.1739e-04, +0.0527, 0.990),
    (4, 7, [3, 5, 1, 6, 0, 2, 4], false, 4.3782e-04, -0.0698, 0.961),
    (3, 5, [1, 3, 0, 4, 2],       false, 4.9499e-04, -0.0993, 0.969),
    (2, 3, [1, 0, 2],             false, 6.0618e-04, -0.0951, 0.907),
    (5, 7, [1, 5, 4, 2, 0, 6, 3], false, 4.0292e-04, -0.0830, 0.981),
    (3, 4, [2, 1, 3, 0],          false, 4.4687e-04, -0.0832, 0.804),
    (4, 5, [2, 3, 1, 4, 0],       true,  3.6578e-04, -0.0494, 0.862),
    (5, 6, [2, 3, 0, 4, 5, 1],    true,  3.2541e-04, -0.0303, 0.861),
    (6, 7, [3, 4, 2, 0, 5, 1, 6], true,  3.0965e-04, -0.0446, 0.989),
    (1, 1, [0],                    false, 2.1100e-04, +0.0029, 1.000),
    (7, 6, [2, 3, 4, 1, 5, 0],    true,  5.7515e-05, +0.1791, 0.970),
    (6, 5, [2, 1, 3, 0, 4],       true,  4.1788e-05, +0.2517, 0.905),
    (5, 4, [2, 1, 0, 3],          true,  6.6505e-05, +0.8915, 0.630),
    (4, 3, [1, 2, 0],             false, 7.8863e-05, +0.5896, 0.909),
    (7, 5, [3, 1, 0, 4, 2],       false, 5.9275e-05, +0.8378, 0.973),
    (3, 2, [0, 1],                false, 1.0672e-04, +0.5522, 0.965),
    (8, 5, [2, 4, 0, 3, 1],       false, 9.0709e-05, +0.9206, 0.988),
    (5, 3, [0, 2, 1],             false, 9.0217e-05, +0.9306, 1.000),
    (7, 4, [0, 3, 2, 1],          false, 8.8271e-05, +0.9573, 0.999),
    (2, 1, [0],                    false, 6.5000e-05, +1.0000, 1.000),
    (7, 3, [2, 0, 1],             false, 1.1053e-05, -0.4062, 0.640),
    (5, 2, [0, 1],                false, 1.5653e-05, -0.2611, 0.733),
    (8, 3, [1, 2, 0],             false, 1.6778e-05, -0.0848, 0.979),
    (3, 1, [0],                    false, 1.7000e-05, +0.0111, 1.000),
]

dnnn_phi = Float64[]
dnnn_abs = Float64[]
dnnn_arg = Float64[]
for dnnn in dNNN_results_80meV
    phi = dnnn[1]/dnnn[2]
    yarg = dnnn[6]
    yabs = dnnn[5]
    if dnnn[end]>0.7 && !(dnnn[4])
        push!(dnnn_phi,phi)
        push!(dnnn_abs,yabs)
        push!(dnnn_arg,yarg)
    end
end


function interpolate_data(data, philist)
    known_data = findall(!ismissing, data)
    order = sortperm(philist[known_data])
    x = Float64.(philist[known_data][order])
    y = Float64.(data[known_data][order])
    return linear_interpolation(x,y)
end

mu_ip = interpolate_data(mu_list, philist)
tnn_ip = interpolate_data(tnn_list, philist)
tnnn_ip = interpolate_data(tnnn_list, philist)
# tdnnn_ip = interpolate_data(tdnnn_list, philist)
tdnnn_abs_ip = interpolate_data(dnnn_abs, dnnn_phi)
tdnnn_arg_ip = interpolate_data(dnnn_arg, dnnn_phi)

tdnnn_arg_ip(1.5)

# kx, ky dimensionless (*a)
function get_h_m(p,q,mu,tnn,tnnn,tdnnn,kx,ky)
    phi = p/q
    h = zeros(ComplexF64, q, q)
    # NN terms
    h += tnn.*diagm([exp(im*(2π*phi*m-ky)) for m = 0:q-1])
    h += tnn.*diagm(-1 => [1 for m = 0:q-2])  
    h[1,q] += tnn*exp(-im*kx)
    # d-NNN terms
    h += tdnnn.*diagm(-1 => [exp(im*(2π*phi*(m+0.5) - ky)) for m = 0:q-2])
    h[1,q] += tdnnn*exp(im*(2π*phi*(q-1+0.5) - kx - ky))
    h += tdnnn.*diagm(1 => [exp(im*(2π*phi*(m-0.5) - ky)) for m = 1:q-1])
    h[q,1] += tdnnn*exp(im*(2π*phi*(0-0.5) + kx - ky))
    # NNN terms
    if q == 1
        h[1,1] += tnnn*exp(-2*im*ky)
        h[1,1] += tnnn*exp(-2*im*kx)
    elseif q == 2
        h += tnnn.*diagm([exp(im*(4π*phi*m-2*ky)) for m = 0:1])
        h += tnnn.*diagm([exp(-im*kx) for m = 0:1])
    elseif q >= 3
        h += tnnn.*diagm([exp(im*(4π*phi*m-2*ky)) for m = 0:q-1])
        h += tnnn.*diagm(-2 => [1 for m = 0:q-3])
        h[2,q] += tnnn*exp(-im*kx)
        h[1,q-1] += tnnn * exp(-im*kx)
    end
    # + h.c.
    h = h + h'
    h .*= -1
    
    # chemical potential
    h += mu.*I
end

p_range = Int(qmax/4):1:2*qmax

# iterate over fluxes
energies = Float64[];
phis = Float64[];
@showprogress for pn in p_range
    # simplify the flux ratio when possible
    gcdpq = gcd(pn,qmax)
    p = div(pn,gcdpq)
    q = div(qmax,gcdpq)
    phi = Float64(p/q)

    kx_range = range(-π/q,π/q,Nk+1)[1:end-1]
    ky_range = range(-π,π,Nk+1)[1:end-1]

    energies_at_phi = Float64[];
    for kx in kx_range
        for ky in ky_range
            h = get_h_m(p,q,mu_ip(phi),tnn_ip(phi),tnnn_ip(phi),tdnnn_abs_ip(phi)*exp(im*π*tdnnn_arg_ip(phi)),kx,ky) 
            append!(energies_at_phi, eigvals(h))
        end
    end
    append!(energies, energies_at_phi)
    append!(phis,[phi for j=1:length(energies_at_phi)])
end



# plot
plt = Plots.scatter(phis,energies,
        xlabel = "", ylabel = "", title = "",
        label = "", ms = 0.8, size = (800,600), color = :yellow,markerstrokewidth=0,
        background_color=:transparent, dpi=200, xaxis=false,yaxis=false, grid=false
);
out_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/modified_butterflies/"
mkpath(out_folder)
savefig(plt,joinpath(out_folder,"modb_v6_U0-80meV.png"))
