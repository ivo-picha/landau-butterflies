module Wannier

using ProgressMeter
    
# function which takes spectrum at given ϕ, finds gaps above critical size, 
# and finds number of states below gap, normalised to particle densities;
# outputs list of particle densities below a gap and a list of corresponding gap sizes
function get_densities_gaps(energy_list::Vector{Float32}, flux::Float32, pn::Int64, NXY::Int64)

    min_gap_size_eV = (energy_list[end] - energy_list[1])/1f2 # minimum gap size to consider significant -- 1% of total plot energy width

    list_big_gaps = Float32[]
    list_densities = Float32[]
    list_energies_gap_lower = Float32[]

    N_en = length(energy_list)
    gaps_list = diff(energy_list)

    for n in eachindex(gaps_list)
        gap = gaps_list[n]
        if gap >= min_gap_size_eV
            push!(list_big_gaps, Float32(gap))
            push!(list_densities, Float32(n))
            push!(list_energies_gap_lower, Float32(energy_list[n])) #such that the upper value of the gap is this + gap
        end
    end
    push!(list_densities, N_en)

    # get some artificial gap size for the last density and gap lower energy
    push!(list_big_gaps, min_gap_size_eV*1f1)
    push!(list_energies_gap_lower, energy_list[end])

    norm_factor = flux/(Float32(pn)*Float32(NXY)^2)
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
const PRECISION = 3 # slope and intercept precision
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
function identify_lines(lines_dict::Dict, unique_phis::Vector{Float32}, datapoints::Vector{NTuple{4, Float32}}, phis_w::Vector{Float32})
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
function merge_round_keys(input_dict::Dict{Tuple{Float32, Float32}, Vector{NTuple{4, Float32}}}; dig::Int = 3)

    merged_dict = Dict{Tuple{Float32, Float32}, Vector{NTuple{4, Float32}}}()

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
function mainW(energies::Vector{Vector{Float32}}, phi_list::Vector{Float32}, qn_list::Vector{Int64}, NXY_list::Vector{Int64})

    num_fluxes = length(phi_list)

    # dictionary to store lines with their points
    lines_dict = Dict{Tuple{Float32, Float32}, Vector{NTuple{4, Float32}}}() # key: (slope, intercept), value: list of points (phi, density, gap size, gap lower energy)

    # loop over fluxes and get densities and gaps at each flux
    all_datapoints = Vector{NTuple{4, Float32}}() # (phi, density, gap size, gap lower energy)
    phis_w = Float32[] # list of phis corresponding to each datapoint

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
            point_tuple = (Float32(phi), Float32(densities_norm[k]), Float32(big_gaps[k]), Float32(energies_gap_lower[k]))
            push!(all_datapoints, point_tuple)
            push!(phis_w, Float32(phi))
        end
    end

    # identify lines in the Wannier plot
    identify_lines(lines_dict, phi_list, all_datapoints, phis_w)

    # merge approximate keys to get integer slope and intercept lines only
    merged_lines_dict = merge_round_keys(lines_dict; dig=PRECISION)

    return merged_lines_dict
    
end
























end