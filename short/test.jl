using Plots
using Measures


spectrum_bare_options = (
    markersize = 0.8,
    color = :black,
    label = "",
    xlabel = "ϕ = p/q",
    ylabel = "E [eV]",
    framestyle = :box,
    size = (1200,800),
    tickfontsize = 16,
    guidefontsize = 18,
    margin = 9mm,
    dpi=900,
)

plt = plot(randn(20);
     spectrum_bare_options...
     )

savefig(plt, "test.png")

