using Plots

# q = [1,2,3,4,5]
# Uc = [0.0151,0.0083,0.00585,0.0047,0.0040]

# scatter([1/qi for qi in q], Uc; 
#     xlabel="1/q", ylabel="Uc", title="Critical potential", markershape=:diamond, label="vLL", color=:red, xlims=(0,1), ylims=(0,0.02))


# scatter(log.(q),log.(Uc); 
#     xlabel="ln q", ylabel="ln Uc", title="Critical potential", markershape=:diamond, label="vLL", color=:red)

# x = log.(q)
# y = log.(Uc)
# X = hcat(ones(length(x)), x)
# coef = X \ y
# intercept, slope = coef[1], coef[2]
# plot!(x, X * coef, label="slope= $(round(slope;digits=4))", color=:red, lw=2)

# plot!(x,X*[coef[1],-1.0])


include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil 
using LinearAlgebra
eigvals(Hamil.get_full_ham(0.1f0, 0f0, 0f0, 0.01f0, 5f-9, 1, 1500))



n(idx,p) = fld(idx, p)
m(idx,p) = idx%p

pp = 5
for idx in 1:(pp*3)
    println("idx=$idx -> n=$(n(idx,pp)), m=$(m(idx,pp))")
end



# visualizing eigenstates

# wf_real = real.(wf)
# wf_imag = imag.(wf)
# wf_abs2 = abs2.(wf)
# wf_arg = angle.(wf)

# wf_real ./= (norm/a^2) # make dimensionless by dividing by a^2
# wf_imag ./= (norm/a^2)
# wf_abs2 ./= (norm/a^2)
# wfr_grid = reshape(wf_real, (Ngrid, Ngrid))
# wfi_grid = reshape(wf_imag, (Ngrid, Ngrid))
# wfa_grid = reshape(wf_abs2, (Ngrid, Ngrid))

## plot wavefunction components
# plt_r = heatmap(xplotrange./a, yplotrange./a, wfr_grid', xlabel="x/a", ylabel="y/a", title="Re(ψ)", aspect_ratio=1, size=(500,500), legend = false)
# plt_i = heatmap(xplotrange./a, yplotrange./a, wfi_grid', xlabel="x/a", ylabel="y/a", title="Im(ψ)", aspect_ratio=1, size=(500,500), legend = false)
# plt_a = heatmap(xplotrange./a, yplotrange./a, sqrt.(wfa_grid)', xlabel="x/a", ylabel="y/a", title="|ψ|", aspect_ratio=1, size=(500,500), legend = false)
# plt_th = heatmap(xplotrange./a, yplotrange./a, reshape(angle.(wf), (Ngrid, Ngrid))', xlabel="x/a", ylabel="y/a", title="arg(ψ)", aspect_ratio=1, size=(500,500), legend = false)
# plt_combo = plot(plt_r, plt_i, plt_a, plt_th, layout = grid(2, 2, hgap = 2mm, vgap = 2mm),
#     aspect_ratio = 1,
#     size = (800,800),
#     margin = 1mm, plot_title="ϕ=$(p)/$(q), U0=$(U0) eV, state #$(ss)")


##
U0 = 0.05
a = 5e-9
e = 1.60217662e-19
m_e = 9.10938356e-31
h = 6.62607015e-34
phi = 1.0

fh = sqrt(4f0 * Float32(π^2) * U0 * e / (m_e * a^2)) # harmonic frequency squared
fc = phi * h / (4f0 * m_e * a^2)

