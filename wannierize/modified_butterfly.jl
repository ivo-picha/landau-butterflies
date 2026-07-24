using Plots
using LinearAlgebra
using Interpolations
using ProgressMeter

# params
Nk = 16 # num k-points
qmax = 120

# fluxes
pqlist = [(0,1),
          (1, 8), (1, 7), (1, 6), (1, 5), (1, 4), (2, 7), (1, 3), (3, 8), (2, 5), (3, 7), 
          (1, 2), (4, 7), (3, 5), (5, 8), (2, 3), (5, 7), (3, 4), (4, 5), (5, 6), (6, 7),
          (7, 8), (1, 1), (8, 7), (7, 6), (6, 5), (5, 4), (4, 3), (7, 5), (3, 2), (8, 5),
          (5, 3), (7, 4), (2, 1), (7, 3), (5, 2), (8, 3), (3, 1)]
philist = [p//q for (p,q) in pqlist]

# data 20 meV
# chemical potential
mu_list = [0.00026486,
       missing, 0.0002836, 0.000284, 0.000304, 0.000313, missing, 0.000359, missing, 0.0004156, missing,
       0.000515, missing, 0.000645, missing, 0.000745, missing, 0.000878, 0.000961, 0.00101, missing,
       missing, 0.001319, missing, 0.00176, missing, missing, 0.002339, missing, 0.003001, 0.003416,
       0.003702, 0.00407, 0.005245, 0.00699, 0.007916, 0.008877, 0.010872]
# NN hopping
tnn_list = [0.00183264,
        missing, 0.001787, 0.00178, 0.001763, 0.001744, missing, 0.001709, missing, 0.001652, missing,
        0.001632, missing, 0.001592, missing, 0.001557, missing, 0.001497, 0.001442, 0.001402, missing,
        missing, 0.001178, missing, 0.000974, missing, missing, 0.000866, missing, 0.000785, 0.000761,
        0.000723, 0.000683, 0.00055, 0.000407, 0.000354, 0.000308, 0.000225]
# NNN hopping
tnnn_list = [0.00022359,
         missing, missing, 0.000205, 0.000203, 0.000203, missing, 0.000201, missing, 0.0002, missing,
         0.00017, missing, 0.00014, missing, 0.000145, missing, 0.000172, 0.000143, 0.000175, missing,
         missing, 0.000123, missing, 6.1e-5, missing, missing, 6.7e-5, missing, 5.8e-5, 5.2e-5,
         4.9e-5, 4.6e-5, 3.3e-5, 1.9e-5, 1.6e-5, 1.2e-5, 7.0e-6] 
# d-NNN hopping
tdnnn_list = [1.4310611595233178e-5,
          missing, missing, 6.0e-5, 9.2e-5, 0.000128, missing, 0.000194, missing, 0.000531, missing,
          0.000323, missing, 0.000454, missing, 0.000979, missing, 0.000382, 0.000318, 0.000241, missing,
          missing, 0.000211, missing, 6.6e-5, missing, missing, 9.1e-5, missing, 0.000167, 0.000106,
          0.000124, 0.0001, 6.5e-5, 1.2e-5, 1.75e-5, 2.2e-5, 1.7e-5]

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
tdnnn_ip = interpolate_data(tdnnn_list, philist)

# kx, ky dimensionless (*a)
function get_h(p,q,mu,tnn,tnnn,tdnnn,kx,ky)
    phi = p/q
    h = zeros(ComplexF64, q, q)
    # NN terms
    h += tnn.*diagm([exp(im*(2π*phi*m-ky)) for m = 0:q-1])
    h += tnn.*diagm(-1 => [exp(-im*kx) for m = 0:q-2])  
    h[1,q] += tnn*exp(-im*kx)
    # d-NNN terms
    h += tdnnn.*diagm(-1 => [exp(im*(2π*phi*(m+0.5) - kx - ky)) for m = 0:q-2])
    h[1,q] += tdnnn*exp(im*(2π*phi*(q-1+0.5) - kx - ky))
    h += tdnnn.*diagm(1 => [exp(im*(2π*phi*(m-0.5) + kx - ky)) for m = 1:q-1])
    h[q,1] += tdnnn*exp(im*(2π*phi*(0-0.5) + kx - ky))
    # NNN terms
    if q == 1
        h[1,1] += tnnn*exp(-2*im*ky)
        h[1,1] += tnnn*exp(-2*im*kx)
    elseif q == 2
        h += tnnn.*diagm([exp(im*(4π*phi*m-2*ky)) for m = 0:1])
        h += tnnn.*diagm([exp(-2*im*kx) for m = 0:1])
    elseif q >= 3
        h += tnnn.*diagm([exp(im*(4π*phi*m-2*ky)) for m = 0:q-1])
        h += tnnn.*diagm(-2 => [exp(-2*im*kx) for m = 0:q-3])
        h[2,q] += tnnn*exp(-2*im*kx)
        h[1,q-1] += tnnn * exp(-2im * kx)
    end
    # + h.c.
    h .*= exp(-im*(phi+0.5)*π)
    h = h + h'
    
    # chemical potential
    h += mu.*I
end

p_range = round(Int64(philist[1])*qmax):1:round(Int64(philist[end])*qmax)

# iterate over fluxes
energies = Float64[];
phis = Float64[];
@showprogress for pn in p_range
    # simplify the flux ratio when possible
    gcdpq = gcd(pn,qmax)
    p = div(pn,gcdpq)
    q = div(qmax,gcdpq)
    phi = Float64(p/q)

    kx_range = range(0,2π/q,Nk+1)[1:end-1]
    ky_range = range(0,2π,Nk+1)[1:end-1]

    energies_at_phi = Float64[];
    for kx in kx_range
        for ky in ky_range
            h = get_h(p,q,mu_ip(phi),tnn_ip(phi),tnnn_ip(phi),tdnnn_ip(phi),kx,ky) 
            append!(energies_at_phi, eigvals(h))
        end
    end
    append!(energies, energies_at_phi)
    append!(phis,[phi for j=1:length(energies_at_phi)])
end



# plot
plt = Plots.scatter(phis,energies,
        xlabel = "ϕ = p/q", ylabel = "Energy [meV]", title = "U₀ = 20 meV",
        label = "", framestyle = :box, ms = 0.5, size = (800,600), color = :black
);
out_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/modified_butterflies/"
mkpath(out_folder)
savefig(plt,joinpath(out_folder,"modb_mexp_U0-20meV.png"))
