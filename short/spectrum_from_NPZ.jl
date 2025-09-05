using NPZ
using Plots
using Measures
using Statistics:mean
using ProgressMeter


folder_path_in = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/mafalda/highres_data_2B/U0.100/"
folder_path_out = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/S_high_res/"

# true: plot Landau levels, false: plot hofstadter bands (bottom q)
plotlandauQ = false;
# 1 or 2; plots respective LL or butterfly band
level = 2;

#specify manually parameters
U0 = 0.100;
a_aa = 50;
q = 240;



read_folder = readdir(folder_path_in)
energies_v = [Vector{Float64}() for _ in 1:length(read_folder)]
phis_v = [Vector{Float64}() for _ in 1:length(read_folder)]

for (fj,file) in enumerate(read_folder)
    filenpz = npzread(joinpath(folder_path_in, file))
    enj = filenpz["y"]
    phij = filenpz["x"]

    phis_v[fj] = phij
    energies_v[fj] = enj
end

energies = reduce(vcat, energies_v)
phis = reduce(vcat, phis_v)

unique_phis = unique(phis)

energies_cut = Float64[];
phis_cut = Float64[];

println("filtering energies...")
@showprogress for phi in unique_phis
    phimask = phis .== phi
    phiens = energies[phimask]
    Nens = length(phiens)
    p = Int(round(phi*q))
    pn = Int(p/gcd(p,q))
    qn = Int(q/gcd(p,q))

    if phi < 1
        Nk = Int(round(Nens/(2qn)))
        if plotlandauQ
            if level == 1
                phiens_cut = phiens[1:Nk*pn]
            elseif level == 2
                phiens_cut = phiens[Nk*pn+1:Nk*2*pn]
            else
                error("level = 1 or 2 only")
            end
        else
            if level == 1
                phiens_cut = phiens[1:Nk*qn]
            elseif level == 2
                phiens_cut = phiens[Nk*qn+1:end]
            else
                error("level = 1 or 2 only")
            end
        end
    else
        Nk = Int(round(Nens/(2pn)))
        if plotlandauQ
            if level == 1
                phiens_cut = phiens[1:Nk*pn]
            elseif level == 2
                phiens_cut = phiens[Nk*pn+1:end]
            else
                error("level = 1 or 2 only")
            end
        else
            if level == 1
                phiens_cut = phiens[1:Nk*qn]
            elseif level == 2
                phiens_cut = phiens[Nk*qn+1:Nk*2*qn]
            else
                error("level = 1 or 2 only")
            end
        end
    end

    append!(energies_cut, phiens_cut)
    append!(phis_cut, [phi for n = 1:length(phiens_cut)])
end

startphi = minimum(phis_cut)
endphi = maximum(phis_cut)

spectrum_bare_options = (
    markersize = 0.25,
    markerstrokewidth = 0.0,
    color = :black,
    label = "",
    xlabel = "ϕ = p/q",
    ylabel = "E [eV]",
    title = plotlandauQ ? 
                        (level == 1 ? "First LL Spectrum; U₀ = $U0 eV, a = $a_aa Å" : "Second LL Spectrum; U₀ = $U0 eV, a = $a_aa Å") :
                        (level == 1 ? "First BB Spectrum; U₀ = $U0 eV, a = $a_aa Å" : "Second BB Spectrum; U₀ = $U0 eV, a = $a_aa Å"),
    framestyle = :box,
    size = (1200,800),
    tickfontsize = 16,
    guidefontsize = 18,
    margin = 9mm
)


plt_s = scatter(phis_cut, energies_cut; spectrum_bare_options...);

#xlims!(plt_s, (-0.02,2.01));
#ylims!(plt_s, (-0.4, 0.0))

#enzf loaded from ZF bands jl
#scatter!(plt_s, zeros(length(enzf)), enzf, ms = 0.5, msw = 0, color = :red, label = "");

plt_name = plotlandauQ ? "highres_npz_sphi$startphi-ephi$endphi-U0$U0-a$a_aa-q$q-LL$level.png" : "highres_npz_sphi$startphi-ephi$endphi-U0$U0-a$a_aa-q$q-BB$level.png"
savefig(plt_s, joinpath(folder_path_out, plt_name))

