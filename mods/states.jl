module States

include(joinpath(@__DIR__,"hamiltonian.jl"))
using .Hamil

using Polynomials, SpecialPolynomials
using LinearAlgebra
using ProgressMeter
using NPZ
using Statistics # used only for one mean()


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
function discard_high_energies_sharp(states_vec_sorted::Vector{Tuple{Float64, Float64, Vector{ComplexF64}}}, phi::Float64, NLL::Int64, np::Float64)
    N_en = length(states_vec_sorted)
    np_max = phi * (NLL+1)
    N_cut = round(Int, N_en * np / np_max)
    if N_cut >= N_en
        println("$np particles per unit cell is more than the allowed maximum ($np_max) at given flux, NLL. Taking $np_max particles per unit cell instead.")
        return states_vec_sorted
    else
        return states_vec_sorted[1:N_cut]
    end
end

# get the el density at xy
function get_el_density_sharp(x::Float64, y::Float64, states_vec_cut::Vector{Tuple{Float64, Float64, Vector{ComplexF64}}}, phi::Float64, a::Float64, p::Int64, NLL::Int64)
    tot = 0.0; # total density at point; to be added

    for state in states_vec_cut
        qnumb_list = [(nj, state[2] + pj*2π/a) for nj=0:NLL for pj=0:(p-1)] # set up the quantum number basis of state
        wf = sum([state[3][i] * landau_lvl_wf(x,y,qnumb_list[i][1],qnumb_list[i][2],phi,a) for i in eachindex(qnumb_list)]) # wavefunction of state at xy
        tot = tot + wf*conj(wf) # add |wf|^2 to total density
    end

    if imag(tot)/real(tot) > 1e-2 && imag(tot) > 1e-5
        println("Imaginary part of the density is non-zero! at position ($x, $y).")
    end
    return Float64(real(tot))
end

# ================= SMEARING / T-dep ===========================

# Fermi-Dirac distribution
function fermi_dirac(en::Float64, eF::Float64, eT::Float64)
    return 1/(exp((en-eF)/eT)+1)
end

# cut off with a bonus so that spectrum can be smeared later without a sharp cutoff
function discard_high_energies_smear(states_vec_sorted::Vector{Tuple{Float64, Float64, Vector{ComplexF64}}}, phi::Float64, NLL::Int64, np::Float64, TeV::Float64)
    N_en = length(states_vec_sorted)
    np_max = phi * (NLL+1)
    N_cut = round(Int, N_en * np / np_max)
    # find position on which energy is suff large
    if N_cut >= N_en
        println("$np particles per unit cell is more than the allowed maximum ($np_max) at given flux, NLL. Taking $np_max particles per unit cell instead. Add more LLs to the calculation!")
        return (states_vec_sorted, states_vec_sorted[end][1])
    else
        EF = states_vec_sorted[N_cut][1] # Fermi energy
        E_cutoff = EF + 4*TeV
        energies_sorted = first.(states_vec_sorted)
        N_cut_plus = searchsortedlast(energies_sorted, E_cutoff; lt = <)
        if N_cut_plus >= N_en
            println("Temperature smear clips the end of the spectrum. Add more LLs to the calculation!")
            return (states_vec_sorted, EF)
        else
            return (states_vec_sorted[1:N_cut_plus], EF) # returns the list of states + fermi energy 
        end
    end
end

# get the el density at xy and smear mid-band states at cutoff
function get_el_density_smear(x::Float64, y::Float64, states_vec_cut::Vector{Tuple{Float64, Float64, Vector{ComplexF64}}}, phi::Float64, a::Float64, p::Int64, NLL::Int64, EF::Float64, TeV::Float64)
    tot = 0.0; # total density at point; to be added

    for state in states_vec_cut
        qnumb_list = [(nj, state[2] + pj*2π/a) for nj=0:NLL for pj=0:(p-1)] # set up the quantum number basis of state
        wf = sum([state[3][i] * landau_lvl_wf(x,y,qnumb_list[i][1],qnumb_list[i][2],phi,a) for i in eachindex(qnumb_list)]) # wavefunction of state at xy
        tot = tot + fermi_dirac(state[1],EF,TeV) * wf*conj(wf) # add |wf|^2 to total density with a weight given by FD distribution
    end

    if imag(tot)/real(tot) > 1e-2 && imag(tot) > 1e-5
        println("Imaginary part of the density is non-zero! at position ($x, $y).")
    end
    return Float64(real(tot))
end
# ============================================================


function get_density_grids(N_uc_x::Int64, N_uc_y::Int64, Nppuc::Int64, states_vec_cut::Vector{Tuple{Float64, Float64, Vector{ComplexF64}}}, phi::Float64, a::Float64, p::Int64, NLL::Int64, np::Float64)
    xgrid = range(a, (N_uc_x+1)*a, N_uc_x*Nppuc)
    ygrid = range(a, (N_uc_y+1)*a, N_uc_y*Nppuc)
    grid_list = reshape(collect(Base.product(xgrid, ygrid)), :)
    density_list = Float64[];
    @showprogress for xy in grid_list
        push!(density_list, get_el_density_sharp(xy[1], xy[2], states_vec_cut, phi, a, p, NLL))
    end
    norm_factor = np/mean(density_list) /a^2
    density_grid = transpose(reshape(norm_factor.*density_list, N_uc_x*Nppuc, N_uc_y*Nppuc))
    return [collect(xgrid), collect(ygrid), Float64.(density_grid)]
end
# add a method to include T --> smearing
function get_density_grids(N_uc_x::Int64, N_uc_y::Int64, Nppuc::Int64, states_vec_cut::Vector{Tuple{Float64, Float64, Vector{ComplexF64}}}, phi::Float64, a::Float64, p::Int64, NLL::Int64, np::Float64, EF::Float64, TeV::Float64)
    xgrid = range(a, (N_uc_x+1)*a, N_uc_x*Nppuc)
    ygrid = range(a, (N_uc_y+1)*a, N_uc_y*Nppuc)
    grid_list = reshape(collect(Base.product(xgrid, ygrid)), :)
    density_list = Float64[];
    @showprogress for xy in grid_list
        push!(density_list, get_el_density_smear(xy[1], xy[2], states_vec_cut, phi, a, p, NLL, EF, TeV))
    end
    norm_factor = np/mean(density_list) /a^2
    density_grid = transpose(reshape(norm_factor.*density_list, N_uc_x*Nppuc, N_uc_y*Nppuc))
    return [collect(xgrid), collect(ygrid), Float64.(density_grid)]
end


    
end