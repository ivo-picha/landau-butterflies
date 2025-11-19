module Aux # contains auxiliary functions and physical constants
export ħ, qe, m_e, kB

using Primes
using ArgParse
    
# physical constants
const ħ = Float32(6.62607015e-34/(2π));  # Planck constant [J s]
const qe = Float32(1.602176634e-19);      # elementary charge [C]
const m_e = Float32(9.1093837139e-31);   # electron mass [kg];
const kB = Float32(8.617333262e-5);       # Boltzmann constant [eV/K]


# interpret parameters from ARGS, with ArgParse.jl ==========================
# container for parsed output
struct ParseResults
    plotQ::Bool
    dataQ::Bool
    wannierQ::Bool
    varyNLLQ::Bool
    plotWQ::Bool
    
    XLL::Int64
    XBF::Int64

    U0::Float32
    a::Float32
    LLmax::Int64
    q::Int64
    phi_s::Float32
    phi_f::Float32
end
# function to parse ARGS and store them in a ParseResults structure
function my_parse(args::Vector{String})
    s = ArgParseSettings()
    @add_arg_table! s begin
        "U0"
            help = "Strength of periodic potential [eV]."
            arg_type = Float32
            default = 0.02
        "a"
            help = "Lattice constant of periodic potential [nm]."
            arg_type = Float32
            default = 5.0
        "LLmax"
            help = "Maximum Landau level index to include (NLL)."
            arg_type = Int
            default = 30
        "q"
            help = "Denominator of flux fraction phi = p/q"
            arg_type = Int
            default = 120
        "phi_s"
            help = "Starting value of flux per unit cell."
            arg_type = Float32
            default = 0.25
        "phi_f"
            help = "Final value of flux per unit cell."
            arg_type = Float32
            default = 1.0
        "--XLL"
            help = "Only output the n-th Landau level (composite bands): specify one integer after the flag, e.g. --XLL 3 gives first 3 LLs."
            arg_type = Int
            default = 0
        "--XBF"
            help = "Only output the n-th Hofstadter butterfly (composite bands): specify one integer after the flag, e.g. --XBF 2 gives first 2 butterflies.                    Note that if both --XLL and --XBF are specified, -XLL is used above flux 1 and -XBF below flux 1."
            arg_type = Int
            default = 0
        "--plot", "-p"
            help = "Output plots of the spectrum and/or Wannier states."
            action = :store_true
        "--data", "-d"
            help = "Output data files of the calculated spectra."
            action = :store_true
        "--wannier", "-w"
            help = "Perform Wannier analysis and output Wannier plots, on whatever "
            action = :store_true
        "--varyNLL", "-l"
            help = "Reduce the number of Landau levels used at higher fluxes, where convergence is faster, to save computation time."
            action = :store_true
        "--plotW", "-z"
            help = "Output Wannier plot figure. Will be set to false if --wannier is not set."
            action = :store_true
    end

    parsed_args = ArgParse.parse_args(args, s)
    
    if parsed_args["plotW"] && !parsed_args["wannier"]
        parsed_args["plotW"] = false
    end

    return ParseResults(
        get(parsed_args, "plot", false),
        get(parsed_args, "data", false),
        get(parsed_args, "wannier", false),
        get(parsed_args, "varyNLL", false),
        get(parsed_args, "plotW", false),
        parsed_args["XLL"],
        parsed_args["XBF"],
        parsed_args["U0"],
        parsed_args["a"],
        parsed_args["LLmax"],
        parsed_args["q"],
        parsed_args["phi_s"],
        parsed_args["phi_f"],
    )
end


function print_startup_message(parsed_struct::ParseResults)
    println("Starting calculation of spectrum with parameters:")
    println("U0 (potential strength) = $(parsed_struct.U0) eV")
    println("a (lattice constant) = $(parsed_struct.a) nm")
    println("Maximum Landau level index = $(parsed_struct.LLmax)")
    println("q (sets resolution) = $(parsed_struct.q)")
    println("phi_s (starting flux) = $(parsed_struct.phi_s)")
    println("phi_f (final flux) = $(parsed_struct.phi_f)\n")
    parsed_struct.XLL == 0 && parsed_struct.XBF == 0 ? println("Outputting full spectrum") :
        parsed_struct.XLL != 0 ? println("Outputting first $(parsed_struct.XLL) Landau levels only") :
        println("Outputting first $(parsed_struct.XBF) Hofstadter butterflies only")
    println("\nOptions selected:")
    println("Plot output: $(parsed_struct.plotQ)")
    println("Data output: $(parsed_struct.dataQ)")
    println("Wannier analysis: $(parsed_struct.wannierQ)")
    println("Vary number of LLs with flux: $(parsed_struct.varyNLLQ)")
    println("Plot Wannier diagram: $(parsed_struct.plotWQ)\n-------------------------------\n")
end




# set up list of p values for given q and flux range
function get_p_list(q::Int64, startphi::Float32, endphi::Float32)
    return unique(Int64.(collect(range(round(q*startphi),round(q*endphi)))))
end

# function to get NXY (number of k-points in each direction) depending on qn so that fluxes with smaller number of bands get more k-points
function get_NXY(qn::Int64, q::Int64, maxNXY::Int64=Int64(127))
    return Int64(round(maxNXY* qn^(-log(q,maxNXY))))
end

# function to get number of LLs to use at given flux, if varying NLL with flux
function get_NLL_at_flux(LLmax::Int64, phi_s::Float32, phi::Float32, LLmin::Int64=Int64(10))
    NLL = LLmin
    Emax = (LLmax + Float32(0.5)) * phi_s # dimensionless max energy at starting flux
    while (NLL + Float32(0.5)) * phi * 0.9 < Emax && NLL < LLmax
        NLL += 1
    end
    return NLL
end



function cut_spectrum(XLL,XBF,num_fluxes,phi_list,qn_list,NLL_list,NXY_list,Elist)
    E_out = Vector{Vector{Float32}}(undef, num_fluxes)
    if XLL != 0 && XBF == 0
        println("Cutting spectrum to specified LLs...\n")
        for j in 1:num_fluxes
            qn = qn_list[j]
            pn = Int64(round(qn*phi_list[j]))
            if pn < (NLL_list[j]+1) # cant request more LLs than available at this flux
                E_out[j] = Elist[j][1:XLL*pn*NXY_list[j]^2]
            end
        end
    elseif XLL == 0 && XBF != 0
        println("Cutting spectrum to specified butterflies...\n")
        for j in 1:num_fluxes
            qn = qn_list[j]
            pn = Int64(round(qn*phi_list[j]))
            E_out[j] = Elist[j][1:Int(min(XBF*qn,pn*(NLL_list[j]+1))*NXY_list[j]^2)]
        end
    elseif XLL != 0 && XBF != 0
        println("Cutting spectrum to specified LLs above phi=1 and butterflies below...\n")
        for j in 1:num_fluxes
            qn = qn_list[j]
            pn = Int64(round(qn*phi_list[j]))
            cutoff = phi<1 ? qn : pn # take qn below phi=1 and pn above
            larger = Int(max(min(XLL,(NLL_list[j]+1)), XBF))
            E_out[j] = Elist[j][1:larger*cutoff*NXY_list[j]^2]
        end
    else
        for j in 1:num_fluxes
            E_out[j] = Elist[j]
        end
    end
    return E_out
end





end # end of module
