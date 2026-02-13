# visualize the uncoupled LLL ground state at ϕ = 1/q 

using Plots

#output folder
outfolder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/jacobi_theta/"
mkpath(outfolder)

# parameters
q = 1
X = 0
Y = X

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
xrange = collect(range(-2.0, 2.0, length=300))
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
size=(475,500), xtickfont = tf, ytickfont = tf, guidefont = lf,
legend=false, #clims=(0.2, 0.4)
)
