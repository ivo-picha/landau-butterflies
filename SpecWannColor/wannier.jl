module Wannier
    
# function which takes spectrum at given ϕ, finds gaps above critical size, 
# and finds number of states below gap, normalised to particle densities;
# outputs list of particle densities below a gap and a list of corresponding gap sizes
function get_densities_gaps(energy_list::Vector{Float64}, min_gap_size_eV::Float64, flux::Float64, NLL::Int64)

    list_big_gaps = Float64[]
    list_densities = Float64[]
    list_energies_gap_lower = Float64[]

    sort!(energy_list)

    N_en = length(energy_list)
    gaps_list = diff(energy_list)

    for n in eachindex(gaps_list)
        gap = gaps_list[n]
        if gap >= min_gap_size_eV
            push!(list_big_gaps, gap)
            push!(list_densities, n)
            push!(list_energies_gap_lower, energy_list[n]) #such that the upper value of the gap is this + gap
        end
    end
    push!(list_densities, N_en-1)

    # get some artificial gap size for the last density and gap lower energy
    push!(list_big_gaps, min_gap_size_eV * 5)
    push!(list_energies_gap_lower, energy_list[end])

    norm_factor = (NLL+1)*flux/(N_en-1)
    list_densities_norm = list_densities.*norm_factor

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
    for n in eachindex(unique_phis)
        if n != 1
            # create list of y values at given x
            index_first = searchsortedfirst(phis_w, unique_phis[n])
            index_last = index_first + searchsortedlast(phis_w[index_first:end], unique_phis[n]) - 1
            points_at_n = datapoints[index_first:index_last]

            for p2 in points_at_n
                # initialize with some absurd values
                nearest_p1_dist = [1e6, 1e6]  # List to store the smallest and second smallest distances
                nearest_p1 = [(0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0)]  # List to store the corresponding points
            
                for p1 in points_at_n_prev
                    d_sq = get_points_distance_sq(p1, p2)
            
                    # update the nearest points based on the current distance
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

end