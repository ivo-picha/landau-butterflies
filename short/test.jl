using LinearAlgebra
using SparseArrays
using ProgressMeter
using Plots

function peierls_hamiltonian(L::Int, φ::Float64, t::Float64 = 1.0)
    N = L * L
    H = spzeros(ComplexF64, N, N)

    function idx(x, y)
        return x + L * y + 1  # Julia uses 1-based indexing
    end

    for x in 0:L-1, y in 0:L-1
        i = idx(x, y)

        # +x hopping (right)
        if x < L - 1
            j = idx(x+1, y)
            phase = -π * φ * y
            H[i, j] = -t * exp(im * phase)
            H[j, i] = -t * exp(-im * phase)  # Hermitian conjugate
        end

        # +y hopping (up)
        if y < L - 1
            j = idx(x, y+1)
            phase = π * φ * x
            H[i, j] = -t * exp(im * phase)
            H[j, i] = -t * exp(-im * phase)
        end
    end

    return Hermitian(Matrix(H))  # Convert to dense matrix for inspection
end


q = 32
plist = 1:q

ens = Float64[];
phis = Float64[];

@showprogress for p in plist
    gcdpq = gcd(p,q)
    pj = Int(p/gcdpq)
    qj = Int(q/gcdpq)

    fl = Float64(pj/qj)

    evals = eigvals(peierls_hamiltonian(2qj, fl, 1.0))

    append!(ens, evals)
    append!(phis, [fl for n in eachindex(evals)])

end


scatter(phis, ens, ms = 0.6, markerstrokewidth = 0, label = "", color = :black)


eigvals(peierls_hamiltonian(4,0.5))