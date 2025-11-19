module Plt
    
include(joinpath(@__DIR__,"hamiltonian.jl"))
using .Hamil

using Plots
using Measures
using ProgressMeter
using Colors, ColorSchemes
using Interpolations
using CairoMakie

spectrum_bare_options = (
    markersize = 0.8,
    color = :black,
    label = "",
    xlabel = "ϕ = p/q",
    ylabel = "E [meV]",
    framestyle = :box,
    size = (1200,800),
    tickfontsize = 16,
    guidefontsize = 18,
    margin = 9mm
)

function plot_bare_spectrum(out_energies::Vector{Vector{Float32}}, phi_list::Vector{Float32})
    plt = Plots.Plot()    # empty plot
    for (j,ens) in enumerate(out_energies)
        xs = [phi_list[j] for i = 1:length(ens)]
        Plots.scatter!(plt, xs, ens; spectrum_bare_options...)
    end
    miny, maxy = minimum(vcat(out_energies...)_), maximum(vcat(out_energies...))
    Plots.ylims!(plt, miny - 0.05f0*(maxy - miny), maxy + 0.05f0*(maxy - miny))
    return plt
end



# add LL guides to spectrum plot
function plot_add_LL_guide!(plot_name::Plots.Plot, startphi::Float32, endphi::Float32, a::Float32, NLL::Int64)
    phi_list_denser = range(startphi, endphi, 300)
    for n = 0:NLL
        y_LLenergies = Hamil.E_LL.(n, phi_list_denser, a)
        Plots.plot!(plot_name, phi_list_denser, y_LLenergies, color = :red, lw = 1, label = "");
    end
end

# add LL guides AND enveloping functions to spectrum plot
function plot_add_LL_guide_env!(plot_name::Plots.Plot, startphi::Float32, endphi::Float32, U0::Float32, a::Float32, NLL::Int64)
    phi_list_denser = range(startphi, endphi, 300)
    xi0_list_denser = sqrt(2π)./ sqrt.(phi_list_denser)
    for n = 0:NLL
        y_LLenergies = Hamil.E_LL.(n, xi0_list_denser, a)
        y_env = 2*U0* Hamil.Ty.(n, n, xi0_list_denser)
        y_env_lower = y_LLenergies .- y_env
        y_env_upper = y_LLenergies .+ y_env
        Plots.plot!(plot_name, phi_list_denser, y_env_lower, fillrange = y_env_upper, color = :red, fillalpha = 0.08, lw = 0, label = "");
        Plots.plot!(plot_name, phi_list_denser, y_LLenergies, color = :red, lw = 1, label = "");
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
    ylims = (0,2)
    #yticks = false
)

# minimum number of points in a line to color it
const min_line_points = 7

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


# plot the full wannier diagram from lines dictionary
function plot_wannier_all(lines_dict::Dict{Tuple{Float32, Float32}, Vector{NTuple{4, Float32}}}, endphi::Float32, NLL::Int64, plots_title::String="Wannier Diagram")

    println("Plotting a Wannier diagram from lines dictionary.")

    # collect all datapoints and gaps from the dictionary
    all_points = Vector{NTuple{4,Float32}}()
    gaps_global = Float32[]
    for (_key, pts) in lines_dict
        for pt in pts
            push!(all_points, pt)
            push!(gaps_global, pt[3]) # gap size
        end
    end

    # fallback if no gaps found
    if isempty(gaps_global)
        gaps_global = [1f0]
    end

    # plot uncolored points
    plot_w = Plots.plot(;wannier_bare_options...)
    Plots.title!(plot_w, plots_title)

    # set a limit of Chern numbers to be colored
    max_colored_chern = clamp(Int(NLL) + 3, 5, 10)

    # point size scaling based on gap sizes
    ggmax = maximum(gaps_global)
    ggmin = minimum(gaps_global)
    msmax = 4
    msmin = 0.8
    rescale_gap(x) = ggmax == ggmin ? (msmin + msmax)/2 : msmin + (x - ggmin)*(msmax - msmin)/(ggmax - ggmin)

    # scatter all raw points (gray)
    for point in all_points
        ms_gap = rescale_gap(point[3])
        Plots.scatter!(plot_w, [point[1]], [point[2]], markerstrokewidth = 0, ms = ms_gap, label = "", color = :gray)
    end

    println("Overlaying colored points.")

    phi_list_denser = range(0f0, endphi, 100)

    # overlay colored lines/points from dictionary
    for (line_key, points) in lines_dict
        if length(points) > min_line_points
            chern = round(line_key[1]; digits = 1)
            if isinteger(chern)
                # guideline
                line_y_vals = line_key[2] .+ line_key[1] .* phi_list_denser
                Plots.plot!(plot_w, phi_list_denser, line_y_vals, color = :gray, alpha = 0.4, lw = 0.5, label = "")
                for point in points
                    ms_gap = rescale_gap(point[3])
                    Plots.scatter!(plot_w, [point[1]], [point[2]], markerstrokewidth = 0, ms = 1.1*ms_gap, label = "", color = gradient_color_plasma(Int(chern), max_colored_chern))
                end
            end
        end
    end

    # build color legend
    legend_values = range(-max_colored_chern, max_colored_chern)
    legend_colors = [gradient_color_plasma(val, max_colored_chern) for val in legend_values]

    legend_plot = Plots.plot(legend=false, xaxis = false, size=(130, 800), xlims=(0, 2), ylims=(-max_colored_chern, max_colored_chern), tickfontsize = 16, guidefontsize = 18, margin=9mm)
    Plots.xlabel!(legend_plot, "C")
    Plots.yticks!(legend_plot, legend_values)
    for (i, val) in enumerate(legend_values)
        rect_color = legend_colors[i]
        y_pos = val - 0.5
        Plots.plot!(legend_plot, Shape([0.1, 2, 2, 0.1], [y_pos, y_pos, y_pos + 1, y_pos + 1]), color=rect_color, label="", lw=0)
    end

    plot_w_l = Plots.plot(legend_plot, plot_w, layout = @layout [a{0.03w} b]);
    return plot_w_l
end


# color the plot of the spectrum
# function color_gaps!(plot_spectrum::Plots.Plot, lines_dict::Dict, unique_phis::Vector{Float32}, NLL::Int64)
#     println("Coloring gaps on the spectrum according to Chern number.")

#     # min distance under which lists of points aren't broken in sublists
#     phi_spacings = diff(unique_phis)
#     phi_thresh = maximum(phi_spacings)
#     # set a limit of Chern numbers to be colored
#     max_colored_chern = clamp(Int(NLL)+3, 5, 10)

#     for (line_key, points) in lines_dict
#         if length(points) > min_line_points
#             chern = round(line_key[1]; digits = 2)
#             if isinteger(chern)
#                 gap_lower_energy = Float32[]
#                 gap_upper_energy = Float32[]
#                 phis_overlay = Float32[]
#                 for point in points
#                     push!(gap_lower_energy, Float32.(point[4]))
#                     upper_energy = point[4] + point[3]
#                     push!(gap_upper_energy, Float32.(upper_energy))
#                     push!(phis_overlay, Float32.(point[1]))
#                 end

#                 # sort lists by phi
#                 perm = sortperm(phis_overlay)
#                 phis_overlay = phis_overlay[perm]
#                 gap_lower_energy = gap_lower_energy[perm]
#                 gap_upper_energy = gap_upper_energy[perm]
    
#                 # SPLIT THE LISTS INTO SUBLISTS WITHOUT LARGE ENERGY JUMPS
#                 gaps = diff(phis_overlay)
#                 split_indices = findall(gaps .> phi_thresh)
#                 split_points = vcat(0, split_indices, length(phis_overlay))
#                 phis_sublists = [phis_overlay[split_points[i]+1:split_points[i+1]] for i in 1:length(split_points)-1]
#                 gap_lower_energy_sublists = [gap_lower_energy[split_points[i]+1:split_points[i+1]] for i in 1:length(split_points)-1]
#                 gap_upper_energy_sublists = [gap_upper_energy[split_points[i]+1:split_points[i+1]] for i in 1:length(split_points)-1]
    
#                 for (phis_sublist, gap_lower_energy_sublist, gap_upper_energy_sublist) in zip(phis_sublists, gap_lower_energy_sublists, gap_upper_energy_sublists)
#                     Plots.plot!(plot_spectrum, phis_sublist, gap_lower_energy_sublist, fillrange = Float32.(gap_upper_energy_sublist), color = gradient_color_plasma(Int(chern), max_colored_chern), fillalpha = 0.5, lw = 0, label = "")
#                 end
#             end
#         end
#     end
# end


# simplify the code for the case of equal spacings of points in flux
function color_gaps_eq!(plot_spectrum::Plots.Plot, lines_dict::Dict, unique_phis::Vector{Float32}, NLL::Int64)
    println("Coloring gaps on the spectrum according to Chern number.")

    # min distance under which lists of points aren't broken in sublists
    phi_spacing = round(unique_phis[2] - unique_phis[1]; digits = 4)
    # set a limit of Chern numbers to be colored
    max_colored_chern = clamp(Int(NLL)+3, 5, 10)

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
    
                # SPLIT THE LISTS INTO SUBLISTS WITHOUT LARGE ENERGY JUMPS
                phigaps = diff(phis_overlay)
                split_indices = findall(phigaps .> (phi_spacing*1.02f0))
                split_points = vcat(0, split_indices, length(phis_overlay))
                phis_sublists = [phis_overlay[(split_points[i]+1):split_points[i+1]] for i in 1:length(split_points)-1]
                gap_lower_energy_sublists = [gap_lower_energy[split_points[i]+1:split_points[i+1]] for i in 1:length(split_points)-1]
                gap_upper_energy_sublists = [gap_upper_energy[split_points[i]+1:split_points[i+1]] for i in 1:length(split_points)-1]
    
                for (phis_sublist, gap_lower_energy_sublist, gap_upper_energy_sublist) in zip(phis_sublists, gap_lower_energy_sublists, gap_upper_energy_sublists)
                    Plots.plot!(plot_spectrum, phis_sublist, gap_lower_energy_sublist, fillrange = Float32.(gap_upper_energy_sublist), color = gradient_color_plasma(Int(chern), max_colored_chern), fillalpha = 0.8, lw = 0, label = "")
                end
            end
        end
    end
end


 


# PLOTTING HEATMAPS OF ELECTRONIC densities

# settings for plot
density_options = (
    xlabel = "x/a",
    ylabel = "y/a",
    framestyle = :box,
    size = (800,820),
    tickfontsize = 18,
    guidefontsize = 20,
    margin = 9mm,
    color = :viridis,
    aspect_ratio = 1,
    titlefontsize = 14,
    titlelocation = :center
    #ylims = (-1,1)
    #yticks = false
)


function plot_density(xgrid::Vector{Float32}, ygrid::Vector{Float32}, density_grid::Matrix{Float32}, a::Float32, nu::Number)
    plot1 = Plots.heatmap(xgrid./a, ygrid./a, density_grid; density_options..., clim = (0,6))#clim = (clamp((nu-1.0), 0.0, nu), nu+1.5) 

    # xsteps = xgrid[1]:a:xgrid[end]
    # xticklist = [string(i,"a") for i in eachindex(xsteps)]
    # my_xticks = (xsteps, xticklist)
    # xticks!(plot1, my_xticks)
    # ysteps = ygrid[1]:a:ygrid[end]
    # yticklist = [string(i,"a") for i in eachindex(ysteps)]
    # my_yticks = (ysteps, yticklist)
    # yticks!(plot1, my_yticks)

    return plot1
end





# function that plots a 3D lines plot of the DOS as a function of flux; used in main_S_DOS.jl
function plot_3D_DOS(DOSs::Vector{Tuple{Vector{Float32},Vector{Float32}}}, phis::Vector{Float32})
    figDOS = Makie.Figure()
    ax = Axis3(figDOS[1, 1], ylabel="ϕ = p/q", xlabel="ϵ [eV]", zlabel="DOS", azimuth = 1.4*π, elevation = 0.2*π)
    for j in eachindex(phis)
        zv = DOSs[j][2]
        xv = DOSs[j][1]
        yv = fill(phis[j], length(xv))
        Makie.lines!(ax, xv, yv, zv, color = :blue, alpha = 0.5)
    end
    return figDOS
end





# function that plots the band structure of the spectrum with Makie; used in main_S_bands.jl
function plot_band_structure(xvec::Vector{Float32}, yvec::Vector{Float32}, band_energies::Vector{Matrix{Float32}}, nplotted::Int=0)

    # figure out how many bands to plot: nbands_plot
    nbands_tot = length(band_energies)
    if nplotted !=0 
        if nplotted > nbands_tot
            error("You wanted to plot $nplotted bands when only $nbands_tot exist.")
        else
            nbands_plot = nplotted
        end
    else
        nbands_plot = nbands_tot
    end

    figBS = CairoMakie.Figure()
    ax = Axis3(figBS[1,1], xlabel = "ky0*a", ylabel = "Y", zlabel = "ϵ [eV]")

    for n in 1:nbands_plot
        CairoMakie.surface!(ax, xvec, yvec, band_energies[n])
    end
    
    return figBS
end










end