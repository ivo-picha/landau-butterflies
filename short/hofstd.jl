using Plots
using LinearAlgebra
using ProgressMeter


q = 240
plist = 0:q

NK = 17
listKx = range(0,2π,NK)[1:end-1]
listKy = range(0,2π,NK)[1:end-1]

function get_H(p::Int, q::Int, Kx, Ky)
    H = zeros(ComplexF64,q,q)
    d(m) = 2*cos(2π*p*m/q + Ky)
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
    for Kx in listKx
        for Ky in listKy
            H = get_H(pj,qj,Kx,Ky)
            eH = eigvals(H)
            append!(enphi, eH)
        end
    end

    append!(elist, enphi)
    append!(philist, [phi for i in eachindex(enphi)])

end

p1 = scatter(philist, elist, markersize = 0.4, color = :black, label = "", xlabel = "ϕ = p/q", ylabel = "ε", framestyle = :box);

savefig(p1, "./b1.png")
