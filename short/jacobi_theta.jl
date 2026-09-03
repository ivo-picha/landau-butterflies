# visualize the uncoupled LLL ground state at ϕ = 1/q 

using Plots
using Measures

#output folder
outfolder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/jacobi_theta/"
mkpath(outfolder)

# parameters
q = 1
X = 3π/2
Y = 3π/2

# ---------------------------------------------

function jacobi_theta(z::Complex, τ::Complex)
    nmax = 100
    sum = 0.0 + 0.0*im
    for n in -nmax:nmax
        sum += exp(π*im*n^2*τ + 2π*im*n*z)
    end
    return sum
end

# a=1 wavefunction
function wf(x::Real, y::Real, X::Real, Y::Real, q::Integer)
    z = x + im*y
    Z = X + im*Y
    tau = im * 1.0 / q
    arg = -(z - Z/(2π)) / q

    wf = jacobi_theta(arg, tau)
    wf *= exp(-π/q * (y-Y/(2π))^2)
    return wf
end

# -----------------------------------------------
xrange = collect(range(0, 1.0*q, length=300))
yrange = xrange
psi_array = Array{Complex}(undef, length(xrange), length(yrange))
for (i,x) in enumerate(xrange)
    for (j,y) in enumerate(yrange)
        psi_array[i,j] = wf(x, y, X, Y, q)
    end
end
# plot
tf = font(18)
lf = font(20)
plt = heatmap(xrange, yrange, abs.(psi_array)', aspect_ratio=1,
    xlabel="x/a", ylabel="y/a", title="|ψ| @ ϕ=1/$q, X=Y=$X", 
    size=(500,500), xtickfont = tf, ytickfont = tf, guidefont = lf,
    legend=:right, framestyle=:box, clims=(0.7, 1.0),
    color = :viridis, margins = 5mm
)

display(plt)

savefig(plt, joinpath(outfolder, "psi_phi1_over_$q-X$X-Y$Y.pdf"))