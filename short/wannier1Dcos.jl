using Plots
using LinearAlgebra
using ProgressMeter


a = 1.
G = 2π/a
V0 = -50.

N = 7; # matrix size /2
Nk = 200; # k-space nsteps

krange = range(-G,G,Nk)


function get_Hmn(k::Float64, G::Float64, V0::Float64, N::Int)
    d(m) = (k+m*G)^2
    H = diagm(0 => [d(m) for m = -N:N], 1 => V0/2 * ones(2N), -1 => V0/2 * ones(2N))
    return Hermitian(H)
end

# # plot band structure
# elist = Float64[];
# klist = Float64[];
# @showprogress for k in krange
#     eks = eigvals(get_Hmn(k, G, V0, N))
#     append!(elist, eks)
#     append!(klist, [k for n = eachindex(eks)])
# end
# specplt = scatter(klist, elist, ms = 2, markerstrokewidth = 0, label = "")

# # find wannier func in lowest band and plot it
# Nx = 200
# xrange = range(-3a,3a,Nx)
# R = 0. # wR center
# wR = zeros(ComplexF64, Nx)
# @showprogress for (j,x) in enumerate(xrange)
#     wRx = 0.
#     for k in klist
#         eks, evecs = eigen(get_Hmn(k, G, V0, N))
#         for (i,n) = enumerate(range(-N,N))
#             wRx += evecs[:,1][i] *exp(im*k*(x-R)) * exp(im*n*G*x)
#         end
#     end
#     wR[j] += wRx
# end
# plot(xrange, real.(wR) .*(1/(Nx*Nk*2N)))



function second_derivative(y::Vector, dx::Float64)
    n = length(y)
    d2y = zeros(n)
    for i in 2:n-1
        d2y[i] = (y[i+1] - 2*y[i] + y[i-1]) / dx^2
    end
    return d2y
end

function get_t_mu(V0::Float64, Nx::Int=500)
    xrange = range(-3a,3a,Nx)
    R1 = 0. # wR centers
    R2 = 1.0
    wR1 = zeros(ComplexF64, Nx)
    wR2 = zeros(ComplexF64, Nx)
    #println("Wannierizing...")
    for (j,x) in enumerate(xrange)
        wRx1 = 0.
        wRx2 = 0.
        for k in klist
            eks, evecs = eigen(get_Hmn(k, G, V0, N))
            for (i,n) = enumerate(range(-N,N))
                wRx1 += evecs[:,1][i] *exp(im*k*(x-R1)) * exp(im*n*G*x)
                wRx2 += evecs[:,1][i] *exp(im*k*(x-R2)) * exp(im*n*G*x)
            end
        end
        wR1[j] += wRx1
        wR2[j] += wRx2
    end
    wR1n = wR1 .*(1/(Nx*Nk*2N))
    wR2n = wR2 .*(1/(Nx*Nk*2N))

    # calculate t and μ

    wR2n_pp = second_derivative(real.(wR2n), minimum(diff(collect(xrange))))
    pot = [V0*cos(x) for x in xrange]

    t_kin = sum([conj(wR1n[i]) * wR2n_pp[i]] for i in eachindex(xrange))
    t_pot = sum([conj(wR1n[i]) * wR2n[i] * pot[i]] for i in eachindex(xrange))
    t = (t_kin + t_pot)/Nx

    mu_kin = sum([conj(wR2n[i]) * wR2n_pp[i]] for i in eachindex(xrange))
    mu_pot = sum([conj(wR2n[i]) * wR2n[i] * pot[i]] for i in eachindex(xrange))
    mu = (mu_kin+mu_pot)/Nx

    return (real(t), real(mu))
end


V0list = range(-101,-1,21)

tlist = Float64[];
mulist = Float64[];
@showprogress for V0i in V0list
    tti, mmi = get_t_mu(V0i)
    push!(tlist, tti[1])
    push!(mulist, mmi[1])
end
tt

scatter(V0list, tlist, label = "", markerstrokewidth = 0, ms = 2, xlabel = "V₀ [a.u.]", ylabel = "t [a.u.]", framestyle = :box)
scatter(V0list, mulist, label = "", markerstrokewidth = 0, ms = 2, xlabel = "V₀ [a.u.]", ylabel = "μ [a.u.]", framestyle = :box)


