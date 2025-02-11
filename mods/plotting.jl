module Plt
    
include(joinpath(@__DIR__,"hamiltonian.jl"))
using .Hamil

using Plots
using Measures
using ProgressMeter
using Colors, ColorSchemes
using Interpolations

spectrum_bare_options = (
    markersize = 0.4,
    color = :black,
    label = "",
    xlabel = "ϕ = p/q",
    ylabel = "E [eV]",
    framestyle = :box,
    size = (1200,800),
    tickfontsize = 16,
    guidefontsize = 18,
    margin = 9mm
)

function plot_spectrum_bare(plot_name::Plots.Plot, phi_list::Vector{Float64}, en_list::Vector{Float64}, plots_title::String)
    println("Plotting bare spectrum.")
    spectrum_bare = scatter(plot_name, phi_list, en_list; spectrum_bare_options...)
    title!(spectrum_bare, plots_title)
    return spectrum_bare
end
# mulitple dispatch in case plot limits are provided
function plot_spectrum_bare(plot_name::Plots.Plot, phi_list::Vector{Float64}, en_list::Vector{Float64}, plots_title::String, x_lims::Tuple{Number, Number}, y_lims::Tuple{Number, Number})
    println("Plotting bare spectrum.")
    spectrum_bare = scatter(plot_name, phi_list, en_list; spectrum_bare_options...,
                            xlims = x_lims, ylims = y_lims)
    title!(spectrum_bare, plots_title)
    return spectrum_bare
end
# multiple dispatch in case plot doesn't already exist
function plot_spectrum_bare(phi_list::Vector{Float64}, en_list::Vector{Float64}, plots_title::String, x_lims::Tuple{Number, Number}, y_lims::Tuple{Number, Number})
    println("Plotting bare spectrum.")
    spectrum_bare = scatter(phi_list, en_list; spectrum_bare_options...,
                            xlims = x_lims, ylims = y_lims)
    title!(spectrum_bare, plots_title)
    return spectrum_bare
end
function plot_spectrum_bare(phi_list::Vector{Float64}, en_list::Vector{Float64}, plots_title::String)
    println("Plotting bare spectrum.")
    spectrum_bare = scatter(phi_list, en_list; spectrum_bare_options...)
    title!(spectrum_bare, plots_title)
    return spectrum_bare
end



# add LL guides to spectrum plot
function plot_add_LL_guide(plot_name::Plots.Plot, startphi::Float64, endphi::Float64, a::Float64, NLL::Int64)
    phi_list_denser = range(startphi, endphi, 300)
    xi0_list_denser = sqrt(2π)./ sqrt.(phi_list_denser)
    for n = 0:NLL
        y_LLenergies = Hamil.E_LL.(n, xi0_list_denser, a)
        plot!(plot_name, phi_list_denser, y_LLenergies, color = :red, lw = 1, label = "");
    end
end

# add LL guides AND enveloping functions to spectrum plot
function plot_add_LL_guide_env(plot_name::Plots.Plot, startphi::Float64, endphi::Float64, U0::Float64, a::Float64, NLL::Int64)
    phi_list_denser = range(startphi, endphi, 300)
    xi0_list_denser = sqrt(2π)./ sqrt.(phi_list_denser)
    for n = 0:NLL
        y_LLenergies = Hamil.E_LL.(n, xi0_list_denser, a)
        y_env = 2*U0* Hamil.Ty.(n, n, xi0_list_denser)
        y_env_lower = y_LLenergies .- y_env
        y_env_upper = y_LLenergies .+ y_env
        plot!(plot_name, phi_list_denser, y_env_lower, fillrange = y_env_upper, color = :red, fillalpha = 0.08, lw = 0, label = "");
        plot!(plot_name, phi_list_denser, y_LLenergies, color = :red, lw = 1, label = "");
    end
end






# plotting wannier plot
wannier_bare_options = (
    xlabel = "ϕ = p/q",
    ylabel = "Particle density",
    framestyle = :box,
    size = (1200,800),
    tickfontsize = 16,
    guidefontsize = 18,
    margin = 9mm,
    ylims = (-1,1)
    #yticks = false
)

# minimum number of points in a line to color it
const min_line_points = 10

# function for custom coloring
function gradient_color_plasma(value::Int, max_value::Int) # thanks chatGPT, not actually plasma anymore
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

function plot_wannier_all(datapoints::Vector{NTuple{4, Float64}}, gaps_global::Vector{Float64}, lines_dict::Dict, endphi::Float64, NLL::Int64, plots_title::String)

    println("Plotting an uncolored Wannier diagram.")
    # plot uncolored points
    plot_w = plot(;wannier_bare_options...)
    title!(plot_w, plots_title)

    # set a limit of Chern numbers to be colored
    max_colored_chern = clamp(Int(NLL)+3, 5, 10)

    # set up a function to scale point sizes between msmin and msmax; rescale if you change png resolution
    ggmax = maximum(gaps_global)
    ggmin = minimum(gaps_global)
    msmax = 4
    msmin = 0.8
    rescale_gap(x) = msmin + (x - ggmin)*(msmax - msmin)/(ggmax - ggmin)

    @showprogress for point in datapoints
        ms_gap = rescale_gap(point[3])
        scatter!(plot_w, [point[1]], [point[2]], markerstrokewidth = 0, ms = ms_gap, label = "", color = :gray)
    end

    println("Overlaying colored points.")

    phi_list_denser = range(0, endphi, 100) # list of phis for plotting guidelines


    @showprogress for (line_key, points) in lines_dict
        if length(points) > min_line_points
            chern = round(line_key[1]; digits = 1)
            if isinteger(chern)
                #plot guideline
                line_y_vals = line_key[2] .+ line_key[1].*phi_list_denser
                plot!(plot_w, phi_list_denser, line_y_vals, color = :gray, alpha = 0.4, lw = 0.5, label = "")
                for point in points
                    ms_gap = rescale_gap(point[3])
                    scatter!(plot_w, [point[1]], [point[2]], markerstrokewidth = 0, ms = 1.1*ms_gap, label = "", color = gradient_color_plasma(Int(chern), max_colored_chern))
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

    plot_w_l = plot(legend_plot, plot_w, layout = @layout [a{0.03w} b]);
    return plot_w_l
end

# color the plot of the spectrum
function color_gaps(plot_spectrum::Plots.Plot, lines_dict::Dict, unique_phis::Vector{Float64}, NLL::Int64)

    # min distance under which lists of points aren't broken in sublists
    phi_spacings = diff(unique_phis)
    phi_thresh = maximum(phi_spacings)
    # set a limit of Chern numbers to be colored
    max_colored_chern = clamp(Int(NLL)+3, 5, 10)

    for (line_key, points) in lines_dict
        if length(points) > min_line_points
            chern = round(line_key[1]; digits = 2)
            if isinteger(chern)
                gap_lower_energy = Float64[]
                gap_upper_energy = Float64[]
                phis_overlay = Float64[]
                for point in points
                    push!(gap_lower_energy, Float64.(point[4]))
                    upper_energy = point[4] + point[3]
                    push!(gap_upper_energy, Float64.(upper_energy))
                    push!(phis_overlay, Float64.(point[1]))
                end
    
                # SPLIT THE LISTS INTO SUBLISTS WITHOUT LARGE ENERGY JUMPS
                gaps = diff(phis_overlay)
                split_indices = findall(gaps .> phi_thresh)
                split_points = vcat(0, split_indices, length(phis_overlay))
                phis_sublists = [phis_overlay[split_points[i]+1:split_points[i+1]] for i in 1:length(split_points)-1]
                gap_lower_energy_sublists = [gap_lower_energy[split_points[i]+1:split_points[i+1]] for i in 1:length(split_points)-1]
                gap_upper_energy_sublists = [gap_upper_energy[split_points[i]+1:split_points[i+1]] for i in 1:length(split_points)-1]
    
                for (phis_sublist, gap_lower_energy_sublist, gap_upper_energy_sublist) in zip(phis_sublists, gap_lower_energy_sublists, gap_upper_energy_sublists)
                    plot!(plot_spectrum, phis_sublist, gap_lower_energy_sublist, fillrange = Float64.(gap_upper_energy_sublist), color = gradient_color_plasma(Int(chern), max_colored_chern), fillalpha = 0.5, lw = 0, label = "")
                end
            end
        end
    end
end

 


# PLOTTING HEATMAPS OF ELECTRONIC densities

# settings for plot
density_options = (
    xlabel = "x",
    ylabel = "y",
    framestyle = :box,
    size = (1200,800),
    tickfontsize = 16,
    guidefontsize = 18,
    margin = 9mm
    #ylims = (-1,1)
    #yticks = false
)


function plot_density(xgrid::Vector{Float64}, ygrid::Vector{Float64}, density_grid::Matrix{Float64}, a::Float64)
    plot1 = heatmap(xgrid, ygrid, density_grid.*(a^2); density_options...)

    xsteps = xgrid[1]:a:xgrid[end]
    xticklist = [string(i,"a") for i in eachindex(xsteps)]
    my_xticks = (xsteps, xticklist)
    xticks!(plot1, my_xticks)
    ysteps = ygrid[1]:a:ygrid[end]
    yticklist = [string(i,"a") for i in eachindex(ysteps)]
    my_yticks = (ysteps, yticklist)
    yticks!(plot1, my_yticks)

    return plot1
end


# DOESNT WORK CURRENTLY; includes interpolation
function plot_density_smooth(xgrid::Vector{Float64}, ygrid::Vector{Float64}, density_grid::Matrix{Float64})

    grid_length = 250 # change for denser grid

    xgrid_denser = range(xgrid[1], xgrid[end], length = grid_length)
    ygrid_denser = range(ygrid[1], ygrid[end], length = grid_length)

    itp = interpolate(density_grid, BSpline(Linear()))
    sizemat = size(density_grid)
    xgrid_denser_itp = range(1, sizemat[2], length = grid_length)
    ygrid_denser_itp = range(1, sizemat[1], length = grid_length)
    density_grid_denser = [itp(yi, xi) for xi in xgrid_denser_itp for yi in ygrid_denser_itp] # itp() needs matrix element index!

    plot1 = heatmap(xgrid_denser_itp, ygrid_denser_itp, density_grid_denser; density_options...)
    return plot1
end





end