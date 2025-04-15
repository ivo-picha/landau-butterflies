module Dens

using LinearAlgebra


# Fermi-Dirac distribution
function fermi_dirac(en::Float64, eF::Float64, eT::Float64)
    return 1/(exp((en-eF)/eT)+1)
end

# cut off with a bonus so that spectrum can be smeared later without a sharp cutoff
function discard_high_energies_smear(states_vec_sorted::Vector{Tuple{Float64, Float64, Float64, Vector{ComplexF64}}}, phi::Float64, NLL::Int64, np::Float64, TeV::Float64)
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



# recursively define a hermite polynomial
function hermite_r(x::Float64, n::Int)::Float64
    if n == 0
        return 1.0
    elseif n == 1
        return 2*x
    else
        return 2*x* hermite_r(x, n-1) - 2*(n-1)* hermite_r(x, n-2)
    end
end

# orthogonal hermite function
function hermite_function(x::Float64, n::Int64)
    An = 1/sqrt(2^n * factorial(n) * sqrt(π))
    return An * hermite_r(x,n) * exp(-x^2 /2)
end

function landau_lvl_wf(x::Float64, y::Float64, n::Int, m::Int, ky0::Float64, Ky::Int, phi::Float64, a::Float64, p::Int)::ComplexF64
    lB = a / sqrt(2π * phi)
    ky = ky0 + 2π*m/a + 2π*p*Ky/a
    xix = x/lB - ky*lB
    return 1/a * hermite_function(xix, n) * exp(im * ky * y)
end


function get_density_list(xyplotlist::Vector{Tuple{Float64,Float64}}, state::Tuple{Float64, Float64, Float64, Vector{ComplexF64}}, nm_list::Vector{Tuple{Int,Int}}, phi::Float64, a::Float64, p::Int, NKy::Int=2)
    dens_list = Float64[];
    for xy in xyplotlist
        x = xy[1]
        y = xy[2]
        wf = sum([state[4][j] * exp(-im*Ky*state[3]) * landau_lvl_wf(x,y,nm_list[j][1],nm_list[j][2],state[2],Int(Ky),phi,a,p) for j in eachindex(nm_list) for Ky = -NKy:NKy])
        push!(dens_list, wf*conj(wf))
    end
    return dens_list
end



end