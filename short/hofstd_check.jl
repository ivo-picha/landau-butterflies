using Plots
using LinearAlgebra
using ProgressMeter
using Base.Threads
nt = nthreads()

t1 = 1.0
t2 = t1

q = 120
plist = 0:2q

NK = 33
listKx = range(0,2π,NK)[1:end-1]
listKy = range(0,2π,NK)[1:end-1]

function get_H1(p::Int, q::Int, Kx, Ky)
    H = zeros(ComplexF64,q,q)
    d(m) = 2*cos(2π*p*2*m/q + Ky)
    H += diagm(0 => [d(m) for m = 0:(q-1)], 1 => ones(q-1), -1 => ones(q-1))
    H[1,q] += exp(im*Kx)
    H[q,1] += exp(-im*Kx)
    return Hermitian(H)
end

function get_H2(p::Int, q::Int, Kx, Ky)
    H = zeros(ComplexF64,q,q)
    d(m) = 2*cos(2π*p*(2*m+1)/q + Ky)
    H += diagm(0 => [d(m) for m = 0:(q-1)], 1 => ones(q-1), -1 => ones(q-1))
    H[1,q] += exp(im*Kx)
    H[q,1] += exp(-im*Kx)
    return Hermitian(H)
end

elist = Float64[];
philist = Float64[];
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
            H1 = t1.*get_H1(pj,qj,Kx,Ky)
            H2 = t2.*get_H2(pj,qj,Kx,Ky)
            eH1 = eigvals(H1)
            eH2 = eigvals(H2)
            append!(buffers[tid], eH1)
            append!(buffers[tid], eH2)
        end
    end

    enphi = reduce(vcat, buffers)

    append!(elist, enphi)
    append!(philist, [phi for i in eachindex(enphi)])

end

p1 = scatter(philist, elist, markersize = 0.4, color = :black, label = "", xlabel = "ϕ = p/q", ylabel = "ε", framestyle = :box);

savefig(p1, "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/hofstdt/checkered.png")
