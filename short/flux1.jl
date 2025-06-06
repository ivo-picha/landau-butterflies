using Plots
using Statistics: mean
using ProgressMeter



nu = 0.1 #filling factor (0 to 1 for LLL @ flux 1)
a = 1. # lattice const

Nxy = 32 # number of steps in each X and Y tilde

xrange = range(0,a,Nxy)
XYlist = reshape(collect(Base.product(xrange, xrange.+0.0001)),:)

energy(X::Float64,Y::Float64,a::Float64) = cos(2π*X/a)+cos(2π*Y/a)

function get_ordered_states(XYlist, a, nu)
    state_list = NTuple{3, Float64}[];
    for XY in XYlist
        push!(state_list, (energy(XY[1], XY[2], a), XY[1], XY[2]))
    end
    sort!(state_list; by=first)
    Nst = length(XYlist)
    return state_list[1:Int(round(Nst*nu))]
end

states = get_ordered_states(XYlist, a, nu)

function jacobi_theta(z::ComplexF64,τ::ComplexF64; ϵ::Float64=10e-10, Nmax::Int=50)::ComplexF64
    s = 1. + im*0.
    for n=1:Nmax
        sn = 2 * exp(im*π*τ*n^2) * cos(2π*n*z)
        s += sn
        if abs(real(sn)) < ϵ && abs(imag(sn)) < ϵ
            #println("reached n=$n")
            break
        end
    end
    return s
end

function phisq(x::Float64,y::Float64,X::Float64,Y::Float64,a::Float64)::ComplexF64
    return exp(-π*(y-Y)^2 /a^2) * jacobi_theta(-1/a * (x-X+im*y-im*Y), 0.0+1.0*im)
end

Nuc = 5
plotres = 128
xrangeplot = range(0,Nuc*a,plotres)
xylistplot = reshape(collect(Base.product(xrangeplot, xrangeplot)),:)

density_list = Float64[];
@showprogress for xy in xylistplot
    dens = 0.;
    for state in states
        ph = phisq(xy[1],xy[2],state[2],state[3],a)
        dens += ph*conj(ph)
    end
    push!(density_list,real(dens)/length(states))
end
density_mat = reshape(density_list, plotres, plotres)
heatmap(xrangeplot,xrangeplot,density_mat, clims=(0,1.2), color = :viridis)