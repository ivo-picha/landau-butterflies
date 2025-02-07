module States

include(joinpath(@__DIR__,"hamiltonian.jl"))
using .Hamil

using Polynomials, SpecialPolynomials
using LinearAlgebra


function landau_lvl_wf(x::Float64, y::Float64, n::Int64, ky::Float64, phi::Float64, a::Float64)::ComplexF64
    lB = a / sqrt(2π * phi)
    xix = x/lB - ky*lB
    Hn = Hermite(n)
    An = 1/sqrt(2^n * factorial(n) * sqrt(π))

    return 1/sqrt(lB*a) * An * Hn(xix) * exp(-xix^2 /2) * exp(im * ky * y)
end

function landau_lvl_wf_c(x::Float64, y::Float64, n::Int64, ky::Float64, phi::Float64, a::Float64)::ComplexF64
    lB = a / sqrt(2π * phi)
    xix = x/lB - ky*lB
    Hn = Hermite(n)
    An = 1/sqrt(2^n * factorial(n) * sqrt(π))

    return 1/sqrt(lB*a) * An * Hn(xix) * exp(-xix^2 /2) * exp(-im * ky * y)
end

# cut off spectrum and eigenvectors to only ones below a given filling np 
function discard_high_energies(energies::Vector{Float64}, vectors::Matrix{ComplexF64}, phi::Float64, NLL::Int64, np::Float64)
    N_en = length(energies)
    np_max = phi * (NLL+1)
    N_cut = round(Int, N_en * np / np_max)
    if N_cut >= N_en
        println("$np particles per unit cell is more than the allowed maximum ($np_max) at given flux, NLL. Taking $np_max particles per unit cell instead.")
        return [energies, vectors]
    else
        return [energies[1:N_cut], vectors[:,1:N_cut]]
    end
end

# get total electronic density at position x,y
function get_el_density(x::Float64, y::Float64, vectors::Vector{Vector{ComplexF64}}, ky_list::Vector{Float64}, phi::Float64, a::Float64, p::Int64, NLL::Int64)
    tot = 0.0; # total density at point; to be added
    n_per_ky = Int(length(vectors)/length(ky_list)) # number of states per ky_star after cutoff

    qnumb_list = [(n, ky) for n=0:NLL for ky=1:p] # ordered list of quantum numbers, as in the eigenvectors

    for k in eachindex(ky_list)
        vs = sum(vectors[(k-1)*n_per_ky + 1 : k*n_per_ky])
        vs_c = conj.(vs)

        sum1 = sum([vs[i]*vs_c[j] * landau_lvl_wf(x,y,qnumb_list[i][1],(ky_list[k]+(qnumb_list[i][2])*2π/a),phi,a) * landau_lvl_wf_c(x,y,qnumb_list[j][1],(ky_list[k]+(qnumb_list[j][2])*2π/a),phi,a) 
                    for i in eachindex(qnumb_list) for j in eachindex(qnumb_list)])
        tot = tot + sum1
    end
    return real(tot)
end

    
end