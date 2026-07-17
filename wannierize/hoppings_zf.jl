# spectrum and states for 2d electron gas in cos potential
using LinearAlgebra
using ProgressMeter



const ħ = 6.62607015e-34/(2π);  # Planck constant [J s]
const e = 1.602176634e-19;      # elementary charge [C]
const m_e = 9.1093837139e-31;   # electron mass [kg];


nu = 1.0; # filling factor
a = 5e-9; # lattice constant
U0 = 0.05; # potential strength in eV
m = m_e; # electron/particle mass

G = 2π/a # recip scat vec
EF = (ħ^2 /e)* 2π*(nu/a^2)/(2*m) 

Nk = 1024; # sqrt of number of momentum states in BZ
NBZ = 4; # number of BZs in each dimension / 2

BZ_centers = reshape(collect(Base.product(G.*(-NBZ:NBZ), G.*(-NBZ:NBZ))), :)
BZ_kpoints = reshape(collect(Base.product(range(-G/2,G/2,Nk), range(-G/2,G/2,Nk))), :)

ϵ = G/1000 # numerical error tolerance

function get_Hk(BZ_kpoint::Tuple{Float64,Float64}, BZ_centers::Vector{Tuple{Float64,Float64}}, m::Float64, ϵ::Float64)
    n = length(BZ_centers)
    Hk = zeros(Float32, n, n)

    # diagonal
    pref = ħ^2 / e
    for i in 1:n
        kx = BZ_kpoint[1] + BZ_centers[i][1]
        ky = BZ_kpoint[2] + BZ_centers[i][2]
        Hk[i,i] = pref * (kx^2 + ky^2) / (2m)
    end

    # off-diagonal
    for i in 1:n, j in i+1:n
        if abs(norm(BZ_centers[i] .- BZ_centers[j]) - G) < ϵ
            Hk[i,j] = U0/2
            Hk[j,i] = U0/2
        end
    end

    return Hermitian(Hk)
end

mu    = 0.0
tnnx  = 0.0 + 0.0im
tnny  = 0.0 + 0.0im
tnnnx = 0.0 + 0.0im
tnnny = 0.0 + 0.0im
tnnnd = 0.0 + 0.0im

@showprogress for BZ_kpoint in BZ_kpoints
    kx, ky = BZ_kpoint

    Hk = get_Hk(BZ_kpoint, BZ_centers, m, ϵ)
    en = minimum(eigvals(Hk))

    # R = (0, 0)
    mu += en

    # R = (a, 0) and (0, a)
    tnnx += en * exp(-im * kx * a)
    tnny += en * exp(-im * ky * a)

    # R = (2a, 0) and (0, 2a)
    tnnnx += en * exp(-im * 2kx * a)
    tnnny += en * exp(-im * 2ky * a)

    # R = (a, a)
    tnnnd += en * exp(-im * (kx + ky) * a)
end

Nktot = length(BZ_kpoints)

mu    /= Nktot
tnnx  /= Nktot
tnny  /= Nktot
tnnnx /= Nktot
tnnny /= Nktot
tnnnd /= Nktot

println("mu    = $mu")
println("tnnx  = $tnnx")
println("tnny  = $tnny")
println("tnnnx = $tnnnx")
println("tnnny = $tnnny")
println("tnnnd = $tnnnd")



