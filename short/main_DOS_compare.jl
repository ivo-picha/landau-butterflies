# spectrum and DOS for 2d electron gas in cos potential with and without field

include(joinpath(dirname(@__DIR__),"mods/params.jl"))
using .Params                       # interpret parameters from ARGS; phys consts; set up lists
include(joinpath(dirname(@__DIR__),"mods/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis                   # wavefunctions and operations for them + density plotting
include(joinpath(dirname(@__DIR__),"mods/plotting.jl"))
using .Plt                          # functions for different kinds of plots      
include(joinpath(dirname(@__DIR__),"mods/densities.jl"))
using .Dens

using LinearAlgebra
using ProgressMeter
using Plots
using Statistics


plot_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/DOS/"
#plot_folder = "/users/ivoga/lh-short/plts"

#number of bins in histograms
Nbins = 50

U0 = 0.03 # potential strenght in eV

# FOR FLUX 1
p = 1
q = 1
NLL = 20
phi = p/q                                   # unit flux per unit cell
xi0 = sqrt(2π / phi)

# FOR FLUX 2
p2 = 2
q2 = 1
NLL2 = 20
phi2 = p2/q2                                  # unit flux per unit cell
xi02 = sqrt(2π / phi2)

a_in_angstr = 50
a = a_in_angstr * 1e-10                     # lattice const in meters
m = 9.1093837139e-31; # electron/particle mass

# ============== nonzero field numerical parameters ============

# get lists of ky0 and Y values to iterate over
Nky = 150;                                   # number of ky* points; independent calculations; variation on scale of U0
ky_list = Params.get_ky_list(a, Nky)
NY = Nky;
Y_list = Params.get_Y_list(NY)

# ============== zero field calc ===================
G = 2π/a # recip scat vec
Nk = 200; # sqrt of number of momentum states in BZ
NBZ = 5; # number of BZs in each dimension / 2
BZ_centers = reshape(collect(Base.product(G.*(-NBZ:NBZ), G.*(-NBZ:NBZ))), :)
BZ_kpoints = reshape(collect(Base.product(range(-G/2,G/2,Nk), range(-G/2,G/2,Nk))), :)
ϵ = G/1000 # numerical error tolerance

function get_Hk(BZ_kpoint::Tuple{Float64,Float64}, BZ_centers::Vector{Tuple{Float64,Float64}}, m::Float64, ϵ::Float64)
    # diagonal elements
    Hk = (ħ^2 /e).*diagm([((BZ_kpoint[1] + kc[1])^2 + (BZ_kpoint[2] + kc[2])^2)/(2*m) for kc in BZ_centers])
    # off-diagonal elements
    for i in eachindex(BZ_centers)
        for j in eachindex(BZ_centers)
            if abs(norm(BZ_centers[i] .- BZ_centers[j]) - G) < ϵ
                Hk[i,j] = U0/2
            end
        end
    end
    return Hermitian(Hk)
end

println("Calculating spectrum and eigenstates for zero field...")

energies_zf = Float64[];
ktots = Float64[];
@showprogress for kpoint in BZ_kpoints
    Hk = get_Hk(kpoint, BZ_centers, m, ϵ)
    evals = eigvals(Hk)
    for i = eachindex(evals)
        append!(energies_zf, evals)
    end
end

sort!(energies_zf)

gaps = diff(energies_zf)
gappos = findfirst(x-> x>U0/2, gaps)
energies_zf[gappos-1]
energies_zf_LB = energies_zf[1:gappos]

# make a histogram output out of intput of energy list, for DOS
function histogram_data(data::Vector{Float64}, N::Int, norm::Float64=1.0)
    # Compute min and max of the data
    min_val = data[1]
    max_val = data[end]
    
    # Compute bin edges and centers
    edges = range(min_val, max_val; length=N+1)
    bin_centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])

    # Initialize frequency vector
    counts = zeros(Int, N)

    # Bin the data manually
    for x in data
        if x == max_val
            # put max value in last bin
            counts[end] += 1
        elseif x == min_val
            counts[1] += 1
        else
            bin_index = searchsortedfirst(edges, x)
            counts[bin_index - 1] += 1
        end
    end

    # Normalize frequencies to sum to np
    total = sum(counts)
    frequencies = counts .* (norm / total)

    return (bin_centers, frequencies)
end

bins_zf, freqs_zf = histogram_data(energies_zf_LB, Nbins)


# ================== nonzero field calculations ===============
println("Calculating nonzero field spectrum...")


energies_f = Float64[];

@showprogress for Y in Y_list
    for ky0 in ky_list

        H = Hamil.get_full_ham(xi0, ky0, Y, U0, a, p, NLL)
        evalsH = eigvals(H)
        append!(energies_f, evalsH)
    end
end
sort!(energies_f)
energies_fs = energies_f[1:(q*(Nky-1)*(NY-1))]

bins_f, freqs_f = histogram_data(energies_fs, Nbins)


# ================== nonzero field 2 calculations ===============
println("Calculating nonzero field spectrum...")


energies_f2 = Float64[];

@showprogress for Y in Y_list
    for ky0 in ky_list

        H = Hamil.get_full_ham(xi02, ky0, Y, U0, a, p2, NLL2)
        evalsH = eigvals(H)
        append!(energies_f2, evalsH)
    end
end
sort!(energies_f2)
energies_fs2 = energies_f2[1:(q2*(Nky-1)*(NY-1))]

bins_f2, freqs_f2 = histogram_data(energies_fs2, Nbins)

# ============== plotting ===================



p1 = Plots.plot(bins_zf,freqs_zf,
    label = "ϕ = 0", framestyle=:box, xlabel="E [eV]", ylabel="DOS", color = :red, title = "U₀ = $U0 eV, a = $a_in_angstr Å")

Plots.plot!(p1,bins_f, freqs_f ,label = "ϕ = $p/$q", color = :blue)

Plots.plot!(p1,bins_f2, freqs_f2 ,label = "ϕ = $p2/$q2", color = :green)

ylims!(p1, (0, maximum([freqs_f;freqs_zf])))
xlims!(p1, (minimum([bins_f2[1],bins_zf[1]]), maximum([bins_f2[end],bins_zf[end]])))


#save plot
plot_name = "DOS_compare_U$U0-f$phi-f2$phi2-N$NLL.png"
savefig(p1, joinpath(plot_folder,plot_name))
