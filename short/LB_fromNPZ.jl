using NPZ
using Plots
using Measures


folder_path_in = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/mafalda/LB_BFs/U0.07/"
folder_path_out = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/S_LB_npz/"

#specify manually parameters
U0 = 0.07;
a_aa = 50;
q = 120;

read_folder = readdir(folder_path_in)
energies_v = [Vector{Float64}() for _ in 1:length(read_folder)]
phis_v = [Vector{Float64}() for _ in 1:length(read_folder)]

for (fj,file) in enumerate(read_folder)
    # if fj == 1 || fj == 4
    #     continue
    # end
    filenpz = npzread(joinpath(folder_path_in, file))
    energies_v[fj] = filenpz["y"]
    phis_v[fj] = filenpz["x"]
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
    title = "Lowest Band Spectrum; U₀ = $U0 eV, a = $a_aa Å",
    framestyle = :box,
    size = (1200,800),
    tickfontsize = 16,
    guidefontsize = 18,
    margin = 9mm
)


plt_s = scatter(phis, energies; spectrum_bare_options...);

plt_name = "LB_npz_sphi$startphi-ephi$endphi-U0$U0-a$a_aa-q$q.png"
savefig(plt_s, joinpath(folder_path_out, plt_name))
