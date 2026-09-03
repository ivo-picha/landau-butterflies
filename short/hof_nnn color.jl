# calculate spectrum for extended an Hofstadter model
# lattice constant a=1
using LinearAlgebra
using Plots
using ProgressMeter
using Colors, ColorSchemes

# path to save plot
out_folder = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/hof_nnn"
mkpath(out_folder)

t1 = 1; # NN hopping
t2 = 0.1; # diagonal NNN hopping
Nkmax = 32; # size of k-grid in each dimension
phi1 = 0; # starting flux
phi2 = 2; # final flux; range defined
q_max = 360; # maximum q value, sets resolution

# bloch hamiltonian; kx ∈ [0,2π/(aq)), ky ∈ [0,2π/a) → magnetic momenta
function get_h(kx::Real, ky::Real, p::Integer, q::Integer)
    phi = p/q
    h = zeros(ComplexF64, q, q)
    # NN terms
    h += -t1.*diagm([exp(im*(2π*phi*m-ky)) for m = 0:q-1])
    h += -t1.*diagm(-1 => [exp(-im*kx) for m = 0:q-2])
    h[1,q] += -t1*exp(-im*kx)
    # NNN terms
    h += -t2.*diagm(-1 => [exp(im*(2π*phi*(m+0.5) - kx - ky)) for m = 0:q-2])
    h[1,q] += -t2*exp(im*(2π*phi*(q-1+0.5) - kx - ky))
    h += -t2.*diagm(1 => [exp(im*(2π*phi*(m-0.5) + kx - ky)) for m = 1:q-1])
    h[q,1] += -t2*exp(im*(2π*phi*(0-0.5) + kx - ky))
    # + h.c.
    h = h + h'
    return Hermitian(h)
end

# p values to iterate over
p_range = 1:1:round(phi2*q_max)

# iterate over fluxes
E_list = Vector{Float64}[];
phis = Float64[];
phi_list = Float64[];
qn_list = Int[];
NXY_list = Int[];
@showprogress for pn in p_range
    # simplify the flux ratio when possible
    gcdpq = gcd(pn,q_max)
    p = div(pn,gcdpq)
    q = div(q_max,gcdpq)
    phi = Float64(p/q)

    Nk = clamp(Int(round(q_max/q))+1,1,Nkmax)
    kx_range = range(0,2π/q,Nk+1)[1:end-1]
    ky_range = range(0,2π,Nk+1)[1:end-1]

    energies_at_phi = Float64[];
    for kx in kx_range
        for ky in ky_range
            h = get_h(kx,ky,p,q)
            append!(energies_at_phi,eigvals(h))
        end
    end
    sort!(energies_at_phi)
    append!(E_list, [energies_at_phi])
    #append!(phis,[phi for j=1:length(energies_at_phi)])
    push!(qn_list,q)
    push!(phi_list,phi)
    push!(NXY_list,Nk)
end



## color ; copied from wannier.jl
function get_densities_gaps(energy_list::Vector{Float64}, flux::Float64, pn::Int64, NXY::Int64)

    min_gap_size_eV = (energy_list[end] - energy_list[1])/1f2 # minimum gap size to consider significant -- 1% of total plot energy width

    list_big_gaps = Float64[]
    list_densities = Float64[]
    list_energies_gap_lower = Float64[]

    N_en = length(energy_list)
    gaps_list = diff(energy_list)

    for n in eachindex(gaps_list)
        gap = gaps_list[n]
        if gap >= min_gap_size_eV
            push!(list_big_gaps, Float64(gap))
            push!(list_densities, Float64(n))
            push!(list_energies_gap_lower, Float64(energy_list[n])) #such that the upper value of the gap is this + gap
        end
    end
    push!(list_densities, N_en)

    # get some artificial gap size for the last density and gap lower energy
    push!(list_big_gaps, min_gap_size_eV*10f1)
    push!(list_energies_gap_lower, energy_list[end])

    norm_factor = flux/(Float64(pn)*Float64(NXY)^2)
    list_densities_norm = list_densities.*norm_factor #round.(list_densities.*norm_factor; digits = 3)

    return [list_densities_norm, list_big_gaps, list_energies_gap_lower]
    
end
# function which identifies points on the same line in Wannier plot and calculates slope + intercept
# updates a dictionary with the coordinates of each unique line and coordinates of each point on the line
# support functions
function get_points_distance_sq(p1::Tuple, p2::Tuple)
    d_sq = (p2[1]-p1[1])^2 + (p2[2]-p1[2])^2
    return d_sq
end
const PRECISION = 4 # slope and intercept precision
function get_slope_intercept_tuple(p1::Tuple, p2::Tuple, precision = PRECISION)
    slope = round( (p2[2]-p1[2])/(p2[1]-p1[1]) ; digits = precision)
    intercept = round( p1[2] - slope*p1[1] ; digits = precision)
    
    if slope == -0. # avoid having different line keys for 0.0 and -0.0
        return (0., intercept)
    elseif intercept == -0.
        return (slope, 0.)
    else
        return (slope, intercept)
    end
end
# main func
function identify_lines(lines_dict::Dict, unique_phis::Vector{Float64}, datapoints::Vector{NTuple{4, Float64}}, phis_w::Vector{Float64})
    # points at previous n for first iteration; gets updated later
    points_at_n_prev = datapoints[1:searchsortedlast(phis_w, unique_phis[1])]
    println("Identifying Wannier diagram lines...")
    @showprogress for n in eachindex(unique_phis)
        if n != 1
            # create list of y values at given x
            index_first = searchsortedfirst(phis_w, unique_phis[n])
            index_last = index_first + searchsortedlast(phis_w[index_first:end], unique_phis[n]) - 1
            points_at_n = datapoints[index_first:index_last]

            for p2 in points_at_n
                # initialize with some absurd values
                nearest_p1_dist = [1f6,1f6]  # List to store the smallest and second smallest distances
                nearest_p1 = [(0f0, 0f0, 0f0, 0f0),(0f0, 0f0, 0f0, 0f0)]  # List to store the corresponding points
            
                for p1 in points_at_n_prev
                    d_sq = get_points_distance_sq(p1, p2)
            
                    # update the nearest 2 points based on the current distance
                    if d_sq < nearest_p1_dist[1]
                        # if d_sq is smaller than the smallest distance
                        nearest_p1_dist[2] = nearest_p1_dist[1]
                        nearest_p1[2] = nearest_p1[1]
                        nearest_p1_dist[1] = d_sq
                        nearest_p1[1] = p1
                    elseif d_sq < nearest_p1_dist[2]
                        # if d_sq is smaller than the second smallest distance but larger than the smallest
                        nearest_p1_dist[2] = d_sq
                        nearest_p1[2] = p1
                    end
                end
            
                # add these points to lines_dict without duplicating points
                for p1 in nearest_p1
                    line_key = get_slope_intercept_tuple(p1, p2) 
            
                    if haskey(lines_dict, line_key)
                        # check if p2 is already in the list to avoid duplicates
                        if !(p2 in lines_dict[line_key])
                            push!(lines_dict[line_key], p2)
                        end
                    else
                        # initialize with both points if line does not exist
                        lines_dict[line_key] = [p1, p2]
                    end
                end
            end
            # update list of points at previous n for the next iteration of n
            points_at_n_prev = points_at_n
        end
    end
end
# merge approximate keys of the lines dictionary and only keep the ones with integer slope and intercept
function merge_round_keys(input_dict::Dict{Tuple{Float64, Float64}, Vector{NTuple{4, Float64}}}; dig::Int = 3)

    merged_dict = Dict{Tuple{Float64, Float64}, Vector{NTuple{4, Float64}}}()

    for (key, vec) in input_dict

        rounded_key = (round(key[1]; digits=dig), round(key[2]; digits=dig))

        if haskey(merged_dict, rounded_key)
            append!(merged_dict[rounded_key], vec)
        elseif isinteger(rounded_key[1]) && isinteger(rounded_key[2])
            merged_dict[rounded_key] = vec
        end
    end

    return merged_dict
end
# main function that performs the whole Wannier analysis and outputs a dictionary of lines with their points
function mainW(energies::Vector{Vector{Float64}}, phi_list::Vector{Float64}, qn_list::Vector{Int64}, NXY_list::Vector{Int64})

    num_fluxes = length(phi_list)

    # dictionary to store lines with their points
    lines_dict = Dict{Tuple{Float64, Float64}, Vector{NTuple{4, Float64}}}() # key: (slope, intercept), value: list of points (phi, density, gap size, gap lower energy)

    # loop over fluxes and get densities and gaps at each flux
    all_datapoints = Vector{NTuple{4, Float64}}() # (phi, density, gap size, gap lower energy)
    phis_w = Float64[] # list of phis corresponding to each datapoint

    println("\nCalculating Wannier data points over flux range...\n")
    @showprogress for j in 1:num_fluxes
        phi = phi_list[j]
        qn = qn_list[j]
        pn = Int64(round(qn*phi))
        NXY = NXY_list[j]
        energies_phi = energies[j]

        wg_output = get_densities_gaps(energies_phi, phi, pn, NXY)
        densities_norm = wg_output[1]
        big_gaps = wg_output[2]
        energies_gap_lower = wg_output[3]

        N_points = length(densities_norm)

        for k in 1:N_points
            point_tuple = (Float64(phi), Float64(densities_norm[k]), Float64(big_gaps[k]), Float64(energies_gap_lower[k]))
            push!(all_datapoints, point_tuple)
            push!(phis_w, Float64(phi))
        end
    end

    # identify lines in the Wannier plot
    identify_lines(lines_dict, phi_list, all_datapoints, phis_w)

    # merge approximate keys to get integer slope and intercept lines only
    merged_lines_dict = merge_round_keys(lines_dict; dig=PRECISION)

    return merged_lines_dict
    
end

wannier_dict = mainW(E_list, phi_list, qn_list, NXY_list)




# plot and save 
function gradient_color_plasma(value::Int, max_value::Int) # thanks chatGPT, not actually plasma anymore
    # Ensure the value is between 0 and max_value
    clamped_value = clamp(value, -max_value, max_value)
    # Normalize the value to a range of 0 to 1
    standardized_value = clamped_value / max_value
    
    # Get the corresponding color from the Plasma colormap
    if standardized_value == 0
        color = colorant"green"
    else
        normalized_value = (1 - standardized_value)/2
        color = ColorSchemes.get(ColorSchemes.balance, normalized_value)
    end
    
    # Return the color as an RGB tuple
    return color
end
function color_gaps_eq2_mainfig(lines_dict::Dict, unique_phis::Vector{Float64})
    println("Coloring gaps on the spectrum according to Chern number.")
    plot_spectrum = Plots.plot()    # empty plot

    # min distance under which lists of points aren't broken in sublists
    phi_spacing = round(unique_phis[2] - unique_phis[1]; digits = 4)
    # set a limit of Chern numbers to be colored
    max_colored_chern = 6 #clamp(Int(NLL)+3, 5, 10)
    min_line_points = 5

    for (line_key, points) in lines_dict
        if length(points) > min_line_points
            chern = round(line_key[1]; digits = 2)
            if isinteger(chern)
                gap_lower_energy = Float32[]
                gap_upper_energy = Float32[]
                phis_overlay = Float32[]
                for point in points
                    push!(gap_lower_energy, Float32.(point[4]))
                    upper_energy = point[4] + point[3]
                    push!(gap_upper_energy, Float32.(upper_energy))
                    push!(phis_overlay, Float32.(point[1]))
                end

                # sort lists by phi
                perm = sortperm(phis_overlay)
                phis_overlay = phis_overlay[perm]
                gap_lower_energy = gap_lower_energy[perm]
                gap_upper_energy = gap_upper_energy[perm]

                en_spacing = (maximum(gap_upper_energy)-minimum(gap_lower_energy))/2f1

                for j in 1:(length(phis_overlay)-1)
                    phi1 = phis_overlay[j]
                    phi2 = phis_overlay[j+1]
                    gl1 = gap_lower_energy[j]
                    gl2 = gap_lower_energy[j+1]
                    gu1 = gap_upper_energy[j]
                    gu2 = gap_upper_energy[j+1]
                    
                    if (phi2-phi1 < 2f0*phi_spacing) && (gu2>gl1) && (gl2<gu1) && (abs(gu2-gu1)<en_spacing) && (abs(gl2-gl1)<en_spacing)
                        Plots.plot!(plot_spectrum, [phi1-0.0001,phi2+0.0001], [gl1,gl2], fillrange = [gu1,gu2],
                        color = gradient_color_plasma(Int(chern), max_colored_chern), fillalpha = 0.7, lw = 0, label = "")
                    end
                end
            end
        end
    end
    return plot_spectrum
end

plot_spectrum = color_gaps_eq2_mainfig(wannier_dict, phi_list)
ylims!(plot_spectrum,-4.6,4.6)


energies = vcat(E_list...)
phis = Float64[];
for (j,phi) in enumerate(phi_list)
    append!(phis,[phi for _ in eachindex(E_list[j])])
end

plt = Plots.scatter!(phis,energies,
        xlabel = "ϕ = p/q", ylabel = "Energy", title = "t = $t1, t'' = $t2",
        label = "", framestyle = :box, ms = 0.5, size = (800,600), color = :black,
        guidefont = 15, labelfont=15, tickfont=13,dpi=2000);
savefig(plt,joinpath(out_folder,"hof_nnn_colored_t1-$t1-t2-$t2-Nk-$Nkmax-s$phi1-f$phi2.png"))