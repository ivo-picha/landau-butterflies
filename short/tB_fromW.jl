# take as input a card, produced by main_t.jl and read off the different hopping parameters
# output a list of eigenenergies
using Plots

input = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/wannier_out/all_hops-ph1-2-U-0.05-a5.0-LL50-NXY48.txt"


# ---------------- read file contents and store values
pots = Vector{Float64}[];
hops = Vector{Float64}[];
p = missing;
q = missing;
U0 = missing;

set = nothing;
open(input, "r") do io
    for (j,line) in enumerate(eachline(io))
        global set
        line = strip(line)
        if j == 2 #read off parameters from second line
            ss = split(line, ", ")
            global U0 = parse(Float64,split(ss[1],"=")[2])
            pq = split(ss[2],"=")[2]
            global p, q = parse.(Int,split(pq,"/"))
            print(ss)
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


# ---------------- build a tight-binding model
function get_tb_h(kx::Real, ky::Real, p::Integer, q::Integer)
    
end