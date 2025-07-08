using Plots
using LinearAlgebra
using ProgressMeter


q = 32
plist = Int.(0:q)

NK = 2
listKx = range(0,2π,NK)[1:end-1]
listKy = range(0,2π,NK)[1:end-1]

function get_Hsymm(p::Int,q::Int,Kx::Float64,Ky::Float64)
    phi = p/q
    nm_list = [(m,n) for m = 1:2q for n = 1:2q]
    H = zeros(ComplexF64,4q^2,4q^2)

    ekx = exp(im*Kx)
    eky = exp(im*Ky)

    # lazy programming
    for (i,nm1) in enumerate(nm_list)
        m1 = nm1[1]
        n1 = nm1[2]
        em = exp(im*π*phi*(m1-1))
        en = exp(im*π*phi*(n1-1))
        for (j,nm2) in enumerate(nm_list)
            m2 = nm2[1]
            n2 = nm2[2]
            if n1 == n2 + 1 && m1 == m2
                H[i,j] += em 
            elseif n1 == n2 - 1 && m1 == m2
                H[i,j] += conj(em) 
            elseif n1 == n2 && m1 == m2 + 1
                H[i,j] += conj(en) 
            elseif n1 == n2 && m1 == m2 - 1
                H[i,j] += en 
            elseif n1 == 2q && n2 == 1 && m1 == m2
                H[i,j] += em * eky
            elseif n1 == 1 && n2 == 2q && m1 == m2
                H[i,j] += conj(em * eky)
            elseif n1 == n2 && m1 == 2q && m2 == 1
                H[i,j] += conj(en) * ekx
            elseif n1 == n2 && m1 == 1 && m2 == 2q
                H[i,j] += en * conj(ekx)
            end
        end
    end
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
            append!(enphi, eigvals(get_Hsymm(pj,qj,Kx,Ky)))
        end
    end
    append!(elist, enphi)
    append!(philist, [phi for i in eachindex(enphi)])
end


p1 = scatter(philist, elist, markersize = 0.4, color = :black, label = "", xlabel = "ϕ = p/q", ylabel = "ε", framestyle = :box);

savefig(p1, "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/hofstdt/symm-NK$(NK-1)-q$q.png")


# krange = range(0,2π,50)
# ens = Float64[];
# for kx in krange
#     for ky in krange
#         en = eigvals(get_Hsymm(1,2,kx,ky))
#         append!(ens, en)
#     end
# end
# histogram(ens, bins=-4:0.1:4)

# round.(get_Hsymm(1,2,0.0,0.0); digits = 4)