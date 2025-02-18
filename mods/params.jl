module Params
using Primes
    
# physical constants
export ħ, e, m_e
const ħ = 6.62607015e-34/(2π);  # Planck constant [J s]
const e = 1.602176634e-19;      # elementary charge [C]
const m_e = 9.1093837139e-31;   # electron mass [kg];

# function to extract parameters from ARGS; called in main_SCWC.jl 
function parse_arguments(args::Vector{String})
    args1 = replace(args[1], "[" => "", "]" => "")
    args2 = split(args1, ",")
    args_vec = parse.(Float64, args2)

    if length(args_vec) != 7
        println("Error: Arguments should be of the format [1.start_field, 2.end_field, 3.U0_in_eV, 4.a_in_Å, 5.p_bands, 6.N_LLs, 7.gap_factor]")
        exit(1)
    end

    if !isprime(Int(args_vec[5]))
        println("\n Code is running but it would be nice if you chose a prime p next time...\n")
    end

    return Any[args_vec[1], args_vec[2], args_vec[3], args_vec[4], Int(args_vec[5]), Int(args_vec[6]), args_vec[7]]
end

# called in main_D.jl file for density calculations
function parse_arguments_D(args::Vector{String})
    args1 = replace(args[1], "[" => "", "]" => "")
    args2 = split(args1, ",")
    args_vec = parse.(Float64, args2)

    if length(args_vec) != 6
        println("Error: Arguments should be of the format [1.p, 2.q, 3.U0_in_eV, 4.a_in_Å, 5.N_LLs, 6.max_particle_density]")
        exit(1)
    end
    return Any[Int(args_vec[1]), Int(args_vec[2]), Float64(args_vec[3]), Float64(args_vec[4]), Int(args_vec[5]), Float64(args_vec[6])]
end

# get a list of q values to interate over when plotting in phi=p/q
function get_q_list(startphi::Number, endphi::Number, p::Int64)
    Nq = Int(round(2*p))
    return unique(Int.(round.([p/(startphi + n*(endphi-startphi)/Nq) for n = range(0, Nq)])))
end

# get a list of q values to interate over when plotting in phi_inv=q/p
function get_q_list_inv(startphi_inv::Number, endphi_inv::Number, p::Int64)
    Nq = Int(round(2*p))
    return unique(Int.(round.(range(startphi_inv * p, endphi_inv * p, Nq))))
end

# get a list of ky* values 
function get_ky_list(a::Float64, Nky::Int64)
    ky_list1 = range(0, 2π/a, Nky)
    return ky_list1[1:end-1]
end

function print_size_message(q_list::Vector, p::Int64, Nky::Int64, NLL::Int64)
    points_per_q = (Nky-1)*(NLL+1)*p
    Nq = length(q_list)
    num_points_in_plot = points_per_q*Nq
    println("Calculating spectrum for $points_per_q x $Nq = $num_points_in_plot points.\n")
end

# extract value of "a" from a filename
function get_a_value(filename::String)
    m = match(r"-a([0-9]+\.[0-9]+)", filename)
    return m !== nothing ? parse(Float64, m.captures[1]) : nothing
end


function startmessage_D(p, phi, U0, a_in_angstr, NLL, np)
    println("\n\n==================
    Calculating the total electronic density of the first $(NLL+1) Landau levels
    in a 2D cos potential with strength U=$U0 eV and lattice constant a=$a_in_angstr Å
    at flux $phi with particle density $np. Resolution of the calculation is p=$p.\n")
end

function startmessage_SCWC(startphi, endphi, U0, a_in_angstr, p, NLL)
    println("\n\n==================
    Calculating the spectrum of the first $(NLL+1) Landau levels
    in a 2D cos potential with strength U=$U0 eV and lattice constant a=$a_in_angstr A
    from flux $startphi to flux $endphi. Resolution of the calculation is p=$p.\n")
end

end
