using NPZ
using Plots
using Measures
using Statistics:mean


folder_path_in = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/mafalda/2B/U0.15/"
folder_path_out = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/S_2B_npz/"

#specify manually parameters
U0 = 0.15;
a_aa = 50;
q = 120;

NLL = 80;

read_folder = readdir(folder_path_in)
energies_v = [Vector{Float64}() for _ in 1:length(read_folder)]
phis_v = [Vector{Float64}() for _ in 1:length(read_folder)]

for (fj,file) in enumerate(read_folder)
    filenpz = npzread(joinpath(folder_path_in, file))
    enj = filenpz["y"]
    phij = filenpz["x"]

    if maximum(phij) < 0.025
        continue        
    end

    phis_v[fj] = phij
    energies_v[fj] = enj
end

energies = reduce(vcat, energies_v)
phis = reduce(vcat, phis_v)

startphi = minimum(phis)
endphi = maximum(phis)

spectrum_bare_options = (
    markersize = 0.4,
    markerstrokewidth = 0.0,
    color = :black,
    label = "",
    xlabel = "ϕ = p/q",
    ylabel = "E [eV]",
    title = "Second Band Spectrum; U₀ = $U0 eV, a = $a_aa Å, Nₗ = $NLL",
    framestyle = :box,
    size = (1200,800),
    tickfontsize = 16,
    guidefontsize = 18,
    margin = 9mm
)


plt_s = scatter(phis, energies; spectrum_bare_options...);

xlims!(plt_s, (-0.02,2.01));
#ylims!(plt_s, (-0.4, 0.0))

#enzf loaded from ZF bands jl
scatter!(plt_s, zeros(length(enzf)), enzf, ms = 0.5, msw = 0, color = :red, label = "");

plt_name = "2B_npz_sphi$startphi-ephi$endphi-U0$U0-a$a_aa-q$q-N$NLL.png"
savefig(plt_s, joinpath(folder_path_out, plt_name))
