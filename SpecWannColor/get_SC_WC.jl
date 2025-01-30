using LinearAlgebra
using Plots
using Primes
using ProgressMeter
using Measures
using Polynomials, SpecialPolynomials
using Serialization
using Colors, ColorSchemes


# GET PARAMETERS FROM ARGUMENTS
args = ARGS
#args = ["[0.05, 1.1, 0.001, 50, 53, 3]"] # for testing

# arguments should be in the format [start_phi, end_phi, U0_in_eV, a_in_Å, p_bands, N_LLs, gap_factor, min_num_points_line]

#get arguments in a vector of floats form
args1 = replace(args[1], "[" => "", "]" => "")
args2 = split(args1, ",")
args_vec = parse.(Float64, args2)

if length(args_vec) != 6
    println("Error: Arguments should be of the format [start_phi_inv, end_phi_inv, U0_in_eV, a_in_Å, p_bands, N_LLs]")
    exit(1)
end

plot_save_folder_path = "~/LH_plots/smooth_plots1"
#data_save_folder_path = "/users/ivoga/single_job_outputs/vecs"

# parameters
a_in_angstr = args_vec[4]
a = args_vec[4]*1e-10 # lattice constant [m]
U0 = args_vec[3] # strength of periodic potential [eV]

# physical constants
ħ = 6.62607015e-34/(2π); # Planck constant [J s]
e = 1.602176634e-19; # elementary charge [C]
m_e = 9.1093837139e-31; # electron mass [kg];

#number of bands to consider 
p_set = Int(args_vec[5])
if !isprime(p_set)
    println("Code is running but it would be nice if you chose a prime p next time...")
end

# list of q values to iterate over
Nq = Int(round(2*p_set)) #number of desired steps of q 
startphi = args_vec[1]
endphi = args_vec[2]
q_list = unique(Int.(round.([p_set/(startphi + n*(endphi-startphi)/Nq) for n = range(0, Nq)])))

# maximum number of LL to consider
NLL = Int(args_vec[6])

#set up list of ky* values to iterate over
Nky = 6;
ky_list = range(0, 2π/a, Nky)
ky_list = ky_list[1:end-1] # 0 = 2π ; remove repetition

# wannier params
gap_factor = 0.4; # the minimum size for a gap to be calculated is gap_factor*ΔE_LL(startphi)
min_line_points = Int(round(p_set/25)); # min number of points in a line to color it







# ================================  FUNCTION DEFINITIONS  ==============================

# define an easier to use laguerre polynomial 
function mylaguerre(α::Number, n::Integer, x::Number)
    list1 = push!(zeros(Integer,n),1)
    lag = Laguerre{α}(list1)
    return lag(x)    
end

# LL energy in eV
E_LL(n::Integer, ξ0::Float64) = (n + 0.5) * ħ^2 / (e * m_e * (ξ0 * a / (2π))^2)

# matrix elements (Θ in overleaf)
Tx(n::Integer, m::Integer, ξ0::Float64) = exp(- ξ0^2 / 4) * (sqrt((2^(n+m)))/(sqrt(factorial(n)) * sqrt(factorial(m))))*
    sum([binomial(n,k)*binomial(m,k) * (1/(2^k)) * factorial(k) * (im*ξ0/2)^(n + m - 2*k) for k = 0:minimum([n,m])])

Ty(n::Integer, m::Integer, ξ0::Float64) = exp(- ξ0^2 / 4) * sqrt(factorial(n))/sqrt(factorial(m)) * (-ξ0/sqrt(2))^(m-n) * mylaguerre(m-n, n, ξ0^2 /2)

# work in basis (n,ky) that goes as (0,ky*), (0,ky*+2π/a), ..., (0,ky*+(p-1)2π/a), (1,ky*),.....
# hamiltonian is composed of A(on diagional) and B matrices that are nonzero on the 3 diagonals

# block matrices on the diagonal -- correspond to (n,k)*(n,k') elements
function matA(n::Integer, ξ0::Float64, ky_star::Float64)
    # get elements from above for LL n
    Txn = Tx(n, n, ξ0)
    Txn_conj = Tx(n, n, -ξ0)
    Tyn = Ty(n, n, ξ0)
    En = E_LL(n, ξ0)

    # diagonal elements
    d(j) = En + U0*(Txn * exp(im * ξ0^2 * a * (ky_star + j * 2π / a)/ (2π)) + Txn_conj * exp(-im * ξ0^2 * a * (ky_star + j * 2π / a)/ (2π)))/2

    # elements on the two second diagonals
    f = U0 * Tyn / 2

    #construct matrix
    mat = diagm(0 => [d(j) for j = 1:p_set], 1 => f * ones(p_set-1), -1 => f * ones(p_set-1))
    mat[1, p_set] = f 
    mat[p_set, 1] = f 

    return mat
end

# matrix elements away from diagonal; LL mixing
function matB(n::Integer, m::Integer, ξ0::Float64, ky_star::Float64)
    # get elements from above for LL n and m
    Txnm = Tx(n, m, ξ0)
    Txmn_conj = Tx(m, n, -ξ0)
    Tynm = Ty(n, m, ξ0)
    Tymn = Ty(m, n, ξ0)

    # diagonal elements
    d(j) = U0*(Txnm * exp(im * ξ0^2 * a * (ky_star + j * 2π / a)/ (2π)) + Txmn_conj * exp(-im * ξ0^2 * a * (ky_star + j * 2π / a)/ (2π)))/2

    # elements on the two second diagonals
    f2 = U0 * Tynm / 2
    f1 = U0 * Tymn / 2

    #construct matrix
    mat = diagm(0 => [d(j) for j = 1:p_set], 1 => f1 * ones(p_set-1), -1 => f2 * ones(p_set-1))
    mat[1, p_set] = f2
    mat[p_set, 1] = f1

    return mat
end

# row n of the hamiltonian matrix
function get_ham_row(n::Integer, ξ0::Float64, ky_star::Float64)
    if n == 0
        arrayB2 = reduce(hcat, [matB(n, m, ξ0, ky_star) for m = n+1:NLL])
        return [matA(n, ξ0, ky_star) arrayB2]
    elseif n == NLL
        arrayB1 = reduce(hcat, [matB(n, m, ξ0, ky_star) for m = 0:n-1])
        return [arrayB1 matA(n, ξ0, ky_star)]
    else
        arrayB1 = reduce(hcat, [matB(n, m, ξ0, ky_star) for m = 0:n-1])
        arrayB2 = reduce(hcat, [matB(n, m, ξ0, ky_star) for m = n+1:NLL])
        return [arrayB1 matA(n, ξ0, ky_star) arrayB2]
    end    
end

# full hamiltonian matrix
function get_full_ham(ξ0::Float64, ky_star::Float64)
    rows_vector = [get_ham_row(n, ξ0, ky_star) for n = 0:NLL]
    return reduce(vcat, rows_vector)    
end



# for wannier diagram
function get_densities_gaps(energy_list, min_gap_size_eV, flux)

    list_big_gaps = []
    list_densities = []
    list_energies_gap_lower = []

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

    # get some artificial gap size for the last density and energy
    push!(list_big_gaps, min_gap_size_eV*5)
    push!(list_energies_gap_lower, energy_list[end])

    norm_factor = (NLL+1)*flux/(N_en-1)
    list_densities_norm = list_densities.*norm_factor

    return [list_densities_norm, list_big_gaps, list_energies_gap_lower]
    
end

# functions for color plotting Wannier
function get_points_distance_sq(p1::Tuple, p2::Tuple)
    d_sq = (p2[1]-p1[1])^2 + (p2[2]-p1[2])^2
    return d_sq
end

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


# function gradient_color(value::Int, max_abs_value::Int) # thanks chatGPT
#     # Ensure max_abs_value is positive
#     max_abs_value = abs(max_abs_value)

#     if value == 0
#         # Return green for zero
#         return colorant"green"
#     elseif value > 0
#         # Map positive values to a gradient from dark red to light red
#         normalized_value = clamp(value / max_abs_value, 0, 1)  # Scale to [0, 1]
#         # Interpolate between dark red and light red
#         return RGB(1 - 0.9 * normalized_value, 0.0, 0.0)  # Red channel increases with `normalized_value`
#     else
#         # Map negative values to a gradient from dark blue to light blue
#         normalized_value = clamp(abs(value) / max_abs_value, 0, 1)  # Scale to [0, 1]
#         # Interpolate between dark blue and light blue
#         return RGB(0.0, 0.0, 1 - 0.9 * normalized_value)  # Blue channel increases with `normalized_value`
#     end
# end

function gradient_color_plasma(value::Number, max_value::Number) # thanks chatGPT, not actually plasma anymore
    # Ensure the value is between 0 and max_value
    clamped_value = clamp(value, -max_value, max_value)
    
    # Normalize the value to a range of 0 to 1
    standardized_value = clamped_value / max_value
    
    # Get the corresponding color from the Plasma colormap
    if standardized_value == 0
        color = colorant"green"
    else
        normalized_value = (1 + standardized_value)/2
        color = get(ColorSchemes.hsv, normalized_value)
    end
    
    # Return the color as an RGB tuple
    return color
end


function merge_approximate_keys(dict::Dict{Tuple{Float64, Float64}, Vector{Tuple{Float64, Float64, Float64, Float64, Float64}}}, tol::Float64)
    # Create a new dictionary to store merged results
    merged_dict = Dict{Tuple{Float64, Float64}, Vector{Tuple{Float64, Float64, Float64, Float64, Float64}}}()
    
    # Track processed keys to avoid duplication
    processed_keys = Set{Tuple{Float64, Float64}}()
    
    for (key1, vec1) in dict
        if key1 in processed_keys
            continue  # Skip already processed keys
        end
        
        # Initialize a merged vector and key for the current group
        merged_vector = vec1
        merged_key = key1
        
        # Find all approximately equal keys
        for (key2, vec2) in dict
            if key2 in processed_keys || key1 == key2
                continue
            end
            
            # Check if keys are approximately equal
            if abs(key1[1] - key2[1]) <= tol && abs(key1[2] - key2[2]) <= tol
                # Merge vectors without duplication
                merged_vector = unique(vcat(merged_vector, vec2))
                # Update the merged key as the mean of the approximately equal keys
                merged_key = ((merged_key[1] + key2[1]) / 2, (merged_key[2] + key2[2]) / 2)
                # Mark key2 as processed
                push!(processed_keys, key2)
            end
        end
        
        # Add the merged result to the new dictionary
        merged_dict[merged_key] = merged_vector
        # Mark key1 as processed
        push!(processed_keys, key1)
    end
    
    return merged_dict
end

# =============================  END OF FUNCTIONS ==================================








# wannier parameters
gap_threshold = gap_factor*(E_LL(1, sqrt(2π / startphi))-E_LL(0, sqrt(2π / startphi))) # minimum size of a gap to be listed (in eV); default -- gap between landau levels at first plotted flux, times factor 2


# size message
sizemat = (NLL+1)*p_set
numberofpoints = sizemat*Nky*Nq
println("Calculating spectrum and Wannier diagram for LLs up to n=$NLL, including $numberofpoints points.")


# empty lists to store coordinates for plot
energies = [];
phis = [];
densities_global = [];
phis_w = []; # separate list to store coordinates for wannier plots
gaps_global = [];
energies_lower_gaps_global = []; # indexing at what energy at a given flux the gap appears

@showprogress for qn in q_list
    p = p_set
    q = qn
    # gcdpq = gcd(p,q)
    # if gcdpq != 1 && p == q # make sure only coprime fractions enter
    #     p = Integer(p / gcdpq)
    #     q = Integer(q / gcdpq)
    # end

    #magnetic field dependent parameters
    phi = p/q # unit flux per unit cell
    xi0 = sqrt(2π / phi)
    lB = xi0 * a / (2π) # magnetic length

    energies_at_phi = [];

    # find spectrum
    for ky in ky_list

        H = get_full_ham(xi0, ky)
        H = Hermitian(H) # declare hermitian to speed up diagonalisation

        evalsH = eigvals(H)
        
        #global energies_at_phi
        energies_at_phi = [energies_at_phi; evalsH] #add eigenvalues to list of energies

    end

    global energies
    energies = [energies; Float64.(energies_at_phi)]

    global phis
    phis = [phis; [phi for j = 1:length(energies_at_phi)]] #add phi values to list of phis

    # get wannier
    densities, gaps, energies_lower_gaps = get_densities_gaps(energies_at_phi, gap_threshold, phi)

    global densities_global
    densities_global = [densities_global; Float64.(densities)]

    global phis_w
    phis_w = [phis_w; [phi for j = 1:length(densities)]]

    global gaps_global
    gaps_global = [gaps_global; Float64.(gaps)]

    global energies_lower_gaps_global
    energies_lower_gaps_global = [energies_lower_gaps_global; Float64.(energies_lower_gaps)]

end







println("Plotting Wannier diagram")

# rescale gap sizes to marker sizes 0.5 to 4; need to be adjusted with png resolution
ggmax = maximum(gaps_global)
ggmin = minimum(gaps_global)
msmax = 4
msmin = 0.8
rescale_gaps(x) = msmin + (x - ggmin)*(msmax - msmin)/(ggmax - ggmin)
marker_sizes = rescale_gaps.(gaps_global)


# adapt color_wannier.jl ===============================================================================
datapoints = collect(zip(phis_w, Float64.(densities_global), marker_sizes, Float64.(gaps_global), Float64.(energies_lower_gaps_global)))
unique_x = unique(phis_w)

# slope and intercet precision
const PRECISION = 4

# find 2 nearest points (w/ smaller x) at each point
# find slope and intercept of line b/w 2 points (key of dict)
# if key exists, add point to dict w/ add_point()
lines_dict = Dict{Tuple{Float64, Float64}, Vector{Tuple{Float64, Float64, Float64, Float64, Float64}}}()

# points at previous n for first iteration; gets updated later
points_at_n_prev = datapoints[1:searchsortedlast(phis_w, unique_x[1])]
@showprogress for n in eachindex(unique_x)
    if n != 1
        # create list of y values at given x
        index_first = searchsortedfirst(phis_w, unique_x[n])
        index_last = index_first + searchsortedlast(phis_w[index_first:end], unique_x[n]) - 1
        points_at_n = datapoints[index_first:index_last]

        for p2 in points_at_n
            # initialize with some absurd values
            nearest_p1_dist = [1e6, 1e6]  # List to store the smallest and second smallest distances
            nearest_p1 = [(0.0, 0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0, 0.0)]  # List to store the corresponding points
        
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
        global points_at_n_prev
        points_at_n_prev = points_at_n
    end
end

# max_abs_chern = maximum(abs.(unique!(chern_list))) # not used atm
max_colored_chern = clamp(Int(NLL)+3, 5, 10)

# get a factor to normalise slope to chern number
testp1 = datapoints[searchsortedlast(phis_w, unique_x[1])]
testp2 = datapoints[searchsortedlast(phis_w, unique_x[end])]
factor_norm_slope = (NLL+1)/(get_slope_intercept_tuple(testp1,testp2)[1])

println("Plotting an uncolored Wannier diagram")
# plot lines with more than 6 points
plot_wannier_lines = plot(title = string("U₀=$U0 eV,  a=$a_in_angstr Å,  Nₗₗ=$NLL"), framestyle = :box, size = (1200,800), tickfontsize = 16, guidefontsize = 18, margin=9mm)
xlims!(plot_wannier_lines, 0, endphi)
ylims!(plot_wannier_lines, -1/factor_norm_slope, 2/factor_norm_slope)
xlabel!(plot_wannier_lines, "Flux ϕ = p/q")
ylabel!(plot_wannier_lines, "Particle density")

# first plot all points and then on top of them the colored ones
@showprogress for point in datapoints
    scatter!(plot_wannier_lines, [point[1]], [point[2]], markerstrokewidth = 0, ms = point[3], label = "", color = :gray)
end

phi_list_denser = range(0, endphi, 200) # list of phis for plotting lines

println("Plotting a colored Wannier diagram")
@showprogress for (line_key, points) in lines_dict
    if length(points) > 1.5*min_line_points
        chern = round(line_key[1]*factor_norm_slope; digits = 1)
        if isinteger(chern)
            # plot line as guide
            line_y_vals = line_key[2] .+ line_key[1].*phi_list_denser
            plot!(plot_wannier_lines, phi_list_denser, line_y_vals, color = :gray, alpha = 0.4, lw = 0.5, label = "")
            for point in points # plot colored points
                scatter!(plot_wannier_lines, [point[1]], [point[2]], markerstrokewidth = 0, ms = 1.1*point[3], label = "", color = gradient_color_plasma(Int(chern), max_colored_chern))
            end
        end
    end
end



# add a gradient legend
legend_values = range(-max_colored_chern, max_colored_chern)
legend_colors = [gradient_color_plasma(val, max_colored_chern) for val in legend_values]

# Create discrete squares for the legend
legend_plot = plot(legend=false, xaxis = false, size=(130, 800), xlims=(0, 2), ylims=(-max_colored_chern, max_colored_chern), tickfontsize = 16, guidefontsize = 18, margin=9mm)
xlabel!(legend_plot, "C")
yticks!(legend_plot, legend_values)
# Plot each value in the legend as a colored square
for (i, val) in enumerate(legend_values)
    rect_color = legend_colors[i]
    y_pos = val - 0.5  # Center the square on the value
    plot!(legend_plot, Shape([0.1, 2, 2, 0.1], [y_pos, y_pos, y_pos + 1, y_pos + 1]), color=rect_color, label="", lw=0)
end

plot_wannier_lines_legended = plot(legend_plot, plot_wannier_lines, layout = @layout [a{0.03w} b]);


wannier_plot_name = string("SgapC_WpsC_outputWzoom_NLL", NLL, "_U0", U0, "_a", a_in_angstr, "_p", p_set, "_from", startphi, "_to", endphi)
wannier_plot_path = string(joinpath(plot_save_folder_path, wannier_plot_name), ".png")
savefig(plot_wannier_lines_legended, wannier_plot_path)



# =================================  PLOTTING SPECTRUM WITH COLORED GAPS =================================================


Npoints = length(energies)
println("Plotting spectrum with $Npoints number of points; coloring gaps appropriately")

plot1 = scatter(phis, energies, title = string("U₀=$U0 eV,  a=$a_in_angstr Å,  Nₗₗ=$NLL"), markersize = 0.4, color = :black, label = "", 
    xlabel = "ϕ = p/q", ylabel = "E [eV]", framestyle = :box, size = (1200,800), tickfontsize = 16, guidefontsize = 18, margin=9mm); # , xlims = (1,1.25), ylims = (200,300)

xlims!(plot1, 0, endphi)

# plot enveloping functions
xi0_list_denser = sqrt(2π)./ sqrt.(phi_list_denser)

function env_upper(n::Integer, xi0::Float64)
    return E_LL(n, xi0) + 2*U0*Ty(n, n, xi0)
end

function env_lower(n::Integer, xi0::Float64)
    return E_LL(n, xi0) - 2*U0*Ty(n, n, xi0)
end

for n = 0:NLL
    #y_env_upper = env_upper.(n, xi0_list_denser)
    #y_env_lower = env_lower.(n, xi0_list_denser)
    y_LLenergies = E_LL.(n, xi0_list_denser)

    #plot!(plot1, phi_list_denser, y_env_lower, fillrange = y_env_upper, color = :red, fillalpha = 0.08, lw = 0, label = "");
    plot!(plot1, phi_list_denser, y_LLenergies, color = :red, lw = 1, label = "");
end


# merge dictionary of lines to avoid any repetitions
merged_dict = merge_approximate_keys(lines_dict, 0.0001)

# COLOR GAPS
phi_spacings = diff(unique_x)
phi_thresh = maximum(phi_spacings)


for (line_key, points) in merged_dict
    if length(points) > min_line_points
        chern = round(line_key[1]*factor_norm_slope; digits = 1)
        if isinteger(chern)
            gap_lower_energy = Float64[]
            gap_upper_energy = Float64[]
            phis_overlay = Float64[]
            for point in points
                push!(gap_lower_energy, Float64.(point[5]))
                upper_energy = point[5] + point[4]
                push!(gap_upper_energy, Float64.(upper_energy))
                push!(phis_overlay, Float64.(point[1]))
            end

            # SPLIT THE LISTS INTO SUBLISTS WITHOUT LARGE ENERGY JUMPS
            gaps = diff(phis_overlay)
            split_indices = findall(gaps .> phi_thresh)
            #split_indices = findall(x -> !(x in phi_spacings), gaps)
            split_points = vcat(0, split_indices, length(phis_overlay))
            phis_sublists = [phis_overlay[split_points[i]+1:split_points[i+1]] for i in 1:length(split_points)-1]
            gap_lower_energy_sublists = [gap_lower_energy[split_points[i]+1:split_points[i+1]] for i in 1:length(split_points)-1]
            gap_upper_energy_sublists = [gap_upper_energy[split_points[i]+1:split_points[i+1]] for i in 1:length(split_points)-1]

            for (phis_sublist, gap_lower_energy_sublist, gap_upper_energy_sublist) in zip(phis_sublists, gap_lower_energy_sublists, gap_upper_energy_sublists)
                plot!(phis_sublist, gap_lower_energy_sublist, fillrange = Float64.(gap_upper_energy_sublist), color = gradient_color_plasma(chern, Float64(max_colored_chern)), fillalpha = 0.5, lw = 0, label = "")
            end
        end
    end
end


spectrum_plot_name = string("SgapC_WpsC_outputS_NLL", NLL, "_U0", U0, "_a", a_in_angstr, "_p", p_set, "_from", startphi, "_to", endphi)
spectrum_plot_path = string(joinpath(plot_save_folder_path, spectrum_plot_name), ".png")
savefig(plot1, spectrum_plot_path)

# save spectrum with static y limits
y_lim_upper = 1.2 * E_LL(NLL, sqrt(2π/endphi))
y_lim_lower = - y_lim_upper/3

ylims!(plot1, y_lim_lower, y_lim_upper)


spectrumZ_plot_name = string("SgapC_WpsC_outputSZ_NLL", NLL, "_U0", U0, "_a", a_in_angstr, "_p", p_set, "_from", startphi, "_to", endphi)
spectrumZ_plot_path = string(joinpath(plot_save_folder_path, spectrumZ_plot_name), ".png")
savefig(plot1, spectrumZ_plot_path)