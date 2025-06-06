using Plots
using LinearAlgebra
using ProgressMeter
using Base.Threads

t = 1.0
tp = 1.0

q = 120
plist = Int.(0:2q)

NK = 65
listKx = range(0,2π,NK)[1:end-1]
listKy = range(0,2π,NK)[1:end-1]

function get_H(p::Int, q::Int, Kx, Ky, t, tp)
    H = zeros(ComplexF64,q,q)
    d(m) = 2*cos(2π*p*m/q + Ky)
    g1(m) = 2*cos(2π*p*(m+0.5)/q + Ky)
    g2(m) = 2*cos(2π*p*(m-0.5)/q + Ky)
    H += diagm(0 => [t*d(m) for m = 0:(q-1)], 1 => [(t+2*tp*g1(m)) for m = 0:(q-2)], -1 => [(t+2*tp*g2(m)) for m = 0:(q-2)])
    H[1,q] += exp(im*Kx)*(t+2*tp*g2(q-1))
    H[q,1] += exp(-im*Kx)*(t+2*tp*g1(q-1))
    return Hermitian(H)
end

elist = Float64[];
philist = Float64[];
nt = nthreads();
@showprogress for p in plist
    gcdpq = gcd(p,q)
    pj = Int(p/gcdpq)
    qj = Int(q/gcdpq)
    phi = pj/qj
    enphi = Float64[];

    buffers = [Vector{Float64}() for _ in 1:nt]

    @threads for Kx in listKx
        tid = threadid()
        for Ky in listKy
            H = get_H(pj,qj,Kx,Ky,t,tp)
            eH = eigvals(H)
            append!(buffers[tid], eH)
        end
    end

    enphi = reduce(vcat, buffers)

    append!(elist, enphi)
    append!(philist, [phi for i in eachindex(enphi)])

end

p1 = scatter(philist, elist, markersize = 0.4, color = :black, label = "", xlabel = "ϕ = p/q", ylabel = "ε", framestyle = :box);
title!(p1, "t=$t, t'=$tp");

savefig(p1, "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/hofstdt/2NN$(NK-1)-q$q-t$t-tp$tp.png")
