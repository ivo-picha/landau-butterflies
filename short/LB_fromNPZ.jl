using NPZ
using Plots
using Measures


filepath = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/local/LB_spectra/LB_S_U0.05-a50.0-q120-phi0.25-phf2.0.npz"
npzfile = npzread(filepath)

phis = npzfile["x"]
energies = npzfile["y"]

defect = minimum(energies)
dmask = energies .!= defect
ennew = energies[dmask]
phnew = phis[dmask]

filepath2 = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/local/LB_spectra/LB_S_U0.05-a50.0-q72-phi2.0-phf3.0.npz"
npzfile2 = npzread(filepath2)

phis2 = npzfile2["x"]
energies2 = npzfile2["y"]

entot = [ennew; energies2]
phtot = [phnew; phis2]

spectrum_bare_options = (
    markersize = 0.4,
    color = :black,
    label = "",
    xlabel = "ϕ = p/q",
    ylabel = "E [eV]",
    title = "Lowest Band; U₀ = 0.05 eV, a = 5 nm",
    framestyle = :box,
    size = (1200,800),
    tickfontsize = 16,
    guidefontsize = 18,
    margin = 9mm
)


plt_s = scatter(phtot, entot; spectrum_bare_options...);

savefig(plt_s, "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/S_LB/SCWC_fixE_vp_LB_S_U0.05-a50.0-qmix-phi0.25-phf3.0.png")


count1 = 0;
for phi in phnew
    if phi == 59/60
        count1 += 1
    end
end
count1/81