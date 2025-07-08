module Params
using Primes
    
# physical constants
export ħ, e, m_e, kB
const ħ = 6.62607015e-34/(2π);  # Planck constant [J s]
const e = 1.602176634e-19;      # elementary charge [C]
const m_e = 9.1093837139e-31;   # electron mass [kg];
const kB = 8.617333262e-5;       # Boltzmann constant [eV/K]

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

    return (args_vec[1], args_vec[2], args_vec[3], args_vec[4], Int(args_vec[5]), Int(args_vec[6]), args_vec[7])
end

# function to extract parameters from ARGS; called in main_SCWC_DOS.jl 
function parse_arguments_DOS(args::Vector{String})
    args1 = replace(args[1], "[" => "", "]" => "")
    args2 = split(args1, ",")
    args_vec = parse.(Float64, args2)

    if length(args_vec) != 7
        println("Error: Arguments should be of the format [1.start_field, 2.end_field, 3.U0_in_eV, 4.a_in_Å, 5.q_fixed, 6.NminLLs, 7.np]")
        exit(1)
    end

    return (args_vec[1], args_vec[2], args_vec[3], args_vec[4], Int(args_vec[5]), Int(args_vec[6]), args_vec[7])
end

# called in main_Dsm.jl file for density calculations
function parse_arguments_D(args::Vector{String})
    args1 = replace(args[1], "[" => "", "]" => "")
    args2 = split(args1, ",")
    args_vec = parse.(Float64, args2)

    if length(args_vec) != 7
        println("Error: Arguments should be of the format [1.p, 2.q, 3.U0_in_eV, 4.a_in_Å, 5.N_LLs, 6.max_particle_density, 7. T]")
        exit(1)
    end
    return (Int(args_vec[1]), Int(args_vec[2]), Float64(args_vec[3]), Float64(args_vec[4]), Int(args_vec[5]), Float64(args_vec[6]), Float64(args_vec[7]))
end

# get a list of q values to interate over when plotting in phi=p/q
function get_q_list(startphi::Number, endphi::Number, p::Int64)
    Nq = Int(round(2*p))
    return unique(Int.(round.([p/(startphi + n*(endphi-startphi)/Nq) for n = range(0, Nq)])))
end

# get a reduced q list that is more equally spaced
function get_q_list_red(startphi::Number, endphi::Number, p::Int64)
    Nq = Int(round(2*p))
    qlist = unique(Int.(round.([p/(startphi + n*(endphi-startphi)/Nq) for n = range(0, Nq)])))
    philist = p./qlist
    phi_d = diff(philist)
    max_phi_d = maximum(phi_d)
    newqlist = [qlist[1]];
    f=2
    c = 0.0;
    for n in eachindex(phi_d)
        if phi_d[n] >= max_phi_d/f
            push!(newqlist, qlist[n+1])
        elseif c < max_phi_d/f
            c += phi_d[n]
            continue
        else
            push!(newqlist, qlist[n+1])
            c = 0.0
        end
    end
    return Int.(newqlist)
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

# get a list of Y values 
function get_Y_list(NY::Int64)
    Ylist1 = range(0, 2π, NY) 
    return Ylist1[1:end-1]
end

function print_size_message(q_list::Vector, p::Int64, Nky::Int64, NY::Int64, NLL::Int64)
    points_per_q = (Nky-1)*(NY-1)*(NLL+1)*p
    Nq = length(q_list)
    num_points_in_plot = points_per_q*Nq
    println("Calculating spectrum for $points_per_q x $Nq = $num_points_in_plot points.\n")
end

# extract value of "a" from a filename
# function get_a_value(filename::String)
#     m = match(r"-a([0-9]+\.[0-9]+)", filename)
#     return m !== nothing ? parse(Float64, m.captures[1]) : nothing
# end

# function extracting parameters from file name
function extract_params(filename::String)
    pattern = r"p(?<p>[^-]+)-q(?<q>[^-]+)-U(?<U0>[^-]+)-a(?<a_in_angstr>[^-]+)-N(?<NLL>[^-]+)-n(?<np>[^-]+)-T(?<TK>[^.]+)"
    m = match(pattern, filename)
    
    if m !== nothing
        return Dict(
            :p => m[:p],
            :q => m[:q],
            :U0 => m[:U0],
            :a_in_angstr => m[:a_in_angstr],
            :NLL => m[:NLL],
            :np => m[:np],
            :TK => m[:TK]
        )
    else
        error("Filename does not match expected pattern.")
    end
end


function startmessage_D(p, phi, U0, a_in_angstr, NLL, np)
    println("\n\n==================
    Calculating the total electronic density of the first $(NLL+1) Landau levels
    in a 2D cos potential with strength U=$U0 eV and lattice constant a=$a_in_angstr Å
    at flux $phi with particle density $np. Resolution of the calculation is p=$p.\n")
end

function startmessage_Dsm(p, phi, U0, a_in_angstr, NLL, np, TK)
    println("\n\n==================
    Calculating the total electronic density of the first $(NLL+1) Landau levels
    in a 2D cos potential with strength U=$U0 eV and lattice constant a=$a_in_angstr Å
    at flux $phi with particle density $np at T=$TK K. Resolution of the calculation is p=$p.\n")
end

function startmessage_SCWC(startphi, endphi, U0, a_in_angstr, p, NLL)
    println("\n\n==================
    Calculating the spectrum of the first $(NLL+1) Landau levels
    in a 2D cos potential with strength U=$U0 eV and lattice constant a=$a_in_angstr A
    from flux $startphi to flux $endphi. Resolution of the calculation is p=$p.\n")
end

function startmessage_S_fixE(startphi, endphi, U0, a_in_angstr, p, Nmin)
    println("\n\n==================
    Calculating the spectrum for a minimum of $(Nmin+1) Landau levels and a maximum of 14
    in a 2D cos potential with strength U=$U0 eV and lattice constant a=$a_in_angstr A
    from flux $startphi to flux $endphi. Resolution of the calculation is p=$p.\n Number of LLs varies at each flux.")
end

end
