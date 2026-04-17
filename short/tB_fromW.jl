# take as input a card, produced by main_t.jl and read off the different hopping parameters
# output a list of eigenenergies
using Plots
using LinearAlgebra
using ProgressMeter
using Measures

input = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_maf/22mar26_t_limited/all_hops-ph3-2-U-0.05-a5.0-LL25-NXY128.txt"
# args = ARGS
# if length(args) != 1
#     println("USAGE: \$0 <all hops txt file>")
#     exit(1)
# end

# input = args[1]

# ---------------- read file contents and store values
pots = Vector{Float64}[];
hops = Vector{Float64}[];
p = missing;
q = missing;
U0 = missing;

set = nothing;
println("Reading from input file...")
open(input, "r") do io
    for (j,line) in enumerate(eachline(io))
        global set
        line = strip(line)
        if j == 2 #read off parameters from second line
            ss = split(line, ", ")
            global U0 = parse(Float64,split(ss[1],"=")[2])
            pq = split(ss[2],"=")[2]
            global p, q = parse.(Int,split(pq,"/"))
        end
        
        # see if data should start being read and what type
        if line == "#chemical-potentials"
            global pots
            set = pots
            continue
        elseif line == "#hopping-amplitudes"
            global hops
            set = hops
            continue
        elseif line == "##chemical-potentials" || line == "##hopping-amplitudes"
            set = nothing
            continue
        end

        if !isnothing(set) && !isempty(line)
            append!(set, [parse.(Float64,split(line))])
        end
    end
end
# ----------------


# ---------------- build a tight-binding model; optional argumen for limiting range of hopping
# basis is m=0 to m=q-1 XY-eigenstates
# can be sped up by switching momenta and orbital summations; not necessary atm
function get_tb_h(X::Real, Y::Real, q::Integer, maxR::Real=10)
    h = zeros(ComplexF64, q, q)
    h += diagm([pots[j][3] for j=1:q])   # on-site potentials
    for hop in hops
        # read off parameters
        Rxfrom = hop[1]
        Ryfrom = hop[2]
        mfrom = Int(hop[3])
        Rxto = hop[4]
        Ryto = hop[5]
        mto = Int(hop[6])
        tabs = hop[7]
        targ = hop[8]
        # hopping vector
        dR = [Rxto-Rxfrom,Ryto-Ryfrom]
        # add terms
        if norm(dR) <= maxR 
            h[mto+1,mfrom+1] += -tabs*exp(im*targ*2π) * exp(im*(X*dR[2] - Y*dR[1])/q)
        end
    end
    return Hermitian(h)
end
# ----------------

# ---------------- plot density of states comparisons
using KernelDensity

# list of Rmax to be compared
maxRlist = [1,sqrt(2),2,3,Inf]

NXY = 128
Xrange = range(0,2π*q,q*NXY+1)[1:end-1]
Yrange = range(0,2π,NXY+1)[1:end-1]
energies = Array{Float64}(undef, q^2*NXY^2,length(maxRlist))
println("Reconsctructing the DOS...")
for (j,maxR) in enumerate(maxRlist)
    ens = Float64[];
    for X in Xrange
        for Y in Yrange
            append!(ens,eigvals(get_tb_h(X,Y,q,maxR)))
        end
    end
    energies[:,j] = ens
end

plt = Plots.plot(framestyle=:box, xlabel="Energy [eV]", ylabel="DOS", yticks=false,
                margins = 5mm, title = "ϕ = $p/$q, U₀ = $U0 eV")
label_list = ["1", "√2", "2", "3", "∞"]
lnwdth = 2;

for j in eachindex(maxRlist)
    ens = energies[:,j]
    eta = (maximum(ens)-minimum(ens))/70
    kd = kde(ens, bandwidth = eta)

    plot!(plt,kd.x,kd.density, lw = lnwdth, label = "ΔRₘₐₓ = $(label_list[j])")
end

# add original DOS; taken from DOS_at_pq.jl
include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil 
phi = Float32(p/q)
a = 5f-9
NLL = Int(round(30*25*U0/phi))
zipped_list = reshape(collect(Iterators.product(Xrange, Yrange)),:)
energieso = Array{Float32}(undef, q, length(zipped_list))
println("Recalculating original spectrum")
@showprogress for i in 1:q*NXY^2
    X = zipped_list[i][1]
    Y = zipped_list[i][2]
    H = Hamil.get_full_ham(phi, Float32(X), Float32(Y), Float32(U0), a, p, NLL)
    energieso[:,i] = eigvals(H)[1:q]
end
energieso = sort(reshape(energieso, :,))
etao = (maximum(energieso)-minimum(energieso))/70
kdo = kde(energieso, bandwidth = etao)
plot!(plt,kdo.x,kdo.density, label = "original", color = :black, style = :dash)


display(plt)
#savefig(plt, joinpath(dirname(input),"DOS_compare_p-$p-q-$q-U0-$U0.png"))

