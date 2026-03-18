# calculate spectrum for extended an Hofstadter model
# lattice constant a=1
using LinearAlgebra
using Plots
using ProgressMeter

# path to save plot
out_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/hof_nnn"
mkpath(out_folder)

t1 = 1; # NN hopping
t2 = 1; # diagonal NNN hopping
Nk = 16; # size of k-grid in each dimension
phi1 = 0; # starting flux
phi2 = 1; # final flux; range defined
q_max = 120; # maximum q value, sets resolution

# bloch hamiltonian; kx ∈ [0,2π/(aq)), ky ∈ [0,2π/a) → magnetic momenta
function get_h(kx::Real, ky::Real, p::Integer, q::Integer)
    phi = p/q
    h = zeros(ComplexF64, q, q)
    # NN terms
    h += -t1.*diagm([exp(im*(2π*phi*m-ky)) for m = 0:q-1])
    h += -t1.*diagm(-1 => [exp(-im*kx) for m = 0:q-2])
    h[1,q] += -t1*exp(-im*kx)
    # NNN terms
    h += -t2.*diagm(-1 => [exp(im*(2π*phi*(m+0.5) - kx - ky)) for m = 0:q-2])
    h[1,q] += -t2*exp(im*(2π*phi*(q-1+0.5) - kx - ky))
    h += -t2.*diagm(1 => [exp(im*(2π*phi*(m-0.5) + kx - ky)) for m = 1:q-1])
    h[q,1] += -t2*exp(im*(2π*phi*(0-0.5) + kx - ky))
    # + h.c.
    h = h + h'
    return Hermitian(h)
end

# p values to iterate over
p_range = round(phi1*q_max):1:round(phi2*q_max)

# iterate over fluxes
energies = Float64[];
phis = Float64[];
@showprogress for pn in p_range
    # simplify the flux ratio when possible
    gcdpq = gcd(pn,q_max)
    p = div(pn,gcdpq)
    q = div(q_max,gcdpq)
    phi = Float64(p/q)

    kx_range = range(0,2π/q,Nk+1)[1:end-1]
    ky_range = range(0,2π,Nk+1)[1:end-1]

    energies_at_phi = Float64[];
    for kx in kx_range
        for ky in ky_range
            h = get_h(kx,ky,p,q)
            append!(energies_at_phi,eigvals(h))
        end
    end
    append!(energies, energies_at_phi)
    append!(phis,[phi for j=1:length(energies_at_phi)])
end

# plot and save 
plt = Plots.scatter(phis,energies,
        xlabel = "ϕ = p/q", ylabel = "Energy", title = "t₁ = $t1, t₂ = $t2",
        label = "", framestyle = :box, ms = 0.4, size = (800,600), color = :black)
savefig(plt,joinpath(out_folder,"hof_nnn_t1-$t1-t2-$t2-Nk-$Nk-s$phi1-f$phi2.png"))