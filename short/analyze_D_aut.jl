# plot a density plot of the electronic density from a pre-saved .npz file

include(joinpath(dirname(@__DIR__), "mods/params.jl"))
using .Params
include(joinpath(dirname(@__DIR__), "mods/densities.jl"))
using .Dens: square_interpolate_d

using NPZ
using Plots
using Peaks
using Interpolations: CubicSplineInterpolation
using LinearAlgebra: diag
using Statistics: mean, std

npz_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/local/dens_grids_p1-q1-U0.01-a50.0-N5-n1.0-T10.0.npz"
plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/npzdens"

a = 50e-10

data = npzread(npz_path)
data_x = data["x"] ./a
data_y = data["y"] ./a
data_z = data["z"]


# interpolate the data for a more precise analysis
data_x_new, data_y_new, data_z_new = Dens.square_interpolate_d(data_x, data_y, data_z)

# check Nuc visually if necessary
heatmap(data_x_new, data_x_new, data_z_new)
# exract the second main diagonal of a matrix
function antidiag(A::Matrix)
    rows, cols = size(A)
    len = min(rows, cols)
    return [A[i, cols - i + 1] for i in 1:len]
end


Nuc = 2
#function get_d_params(data_x::Vector{Float64}, data_y::Vector{Float64}, data_z::Matrix{Float64}, Nuc = 2)

fm = findmax(data_z_new)
translation_index = round(Int, length(data_x_new)/Nuc)

# set the number of steps for an interpolation of the cross-sections
Nitp = 512

# vertical cross-sections
x_peaks_data = [];
x_interp = range(data_x_new[1], data_x_new[end], Nitp)
for i = 1:Nuc
    cs = data_z_new[Int(mod( fm[2][1] + translation_index*i, translation_index*Nuc)),:]
    itp_i = Interpolations.CubicSplineInterpolation(range(data_x_new[1], data_x_new[end], length(data_x_new)), cs)
    itp_data = itp_i.(x_interp)

    #Peaks.jl for analysing cross-sections
    ind_i, heights_i = findmaxima(itp_data)
    ind_i, proms_i = peakproms(ind_i, itp_data)
    ind_i, widths_i, edges_i = peakwidths(ind_i, itp_data, proms_i)

    widths_i_a = widths_i .* ((data_x_new[end]-data_x_new[1])/length(itp_data))
    for j in eachindex(ind_i)
        push!(x_peaks_data, (proms_i[j], widths_i_a[j]))
    end
end

# vertical cross-sections
y_peaks_data = [];
y_interp = range(data_y_new[1], data_y_new[end], Nitp)
for i = 1:Nuc
    cs = data_z_new[:, Int(mod( fm[2][2] + translation_index*i, translation_index*Nuc))]
    itp_i = Interpolations.CubicSplineInterpolation(range(data_y_new[1], data_y_new[end], length(data_y_new)), cs)
    itp_data = itp_i.(y_interp)

    #Peaks.jl for analysing cross-sections
    ind_i, heights_i = findmaxima(itp_data)
    ind_i, proms_i = peakproms(ind_i, itp_data)
    ind_i, widths_i, edges_i = peakwidths(ind_i, itp_data, proms_i)

    widths_i_a = widths_i .* ((data_y_new[end]-data_y_new[1])/length(itp_data))
    for j in eachindex(ind_i)
        push!(y_peaks_data, (proms_i[j], widths_i_a[j]))
    end
end


# diagonal cross-sections on the two main diagonals 
d_peaks_data = [];
d_interp = range(sqrt(data_x_new[1]^2 + data_y_new[1]^2), sqrt(data_x_new[end]^2 + data_y_new[end]^2), Nitp)
d1 = diag(data_z_new)
d2 = antidiag(data_z_new)
itp_d1 = Interpolations.CubicSplineInterpolation(range(d_interp[1], d_interp[end], length(data_x_new)), d1)
itp_d2 = Interpolations.CubicSplineInterpolation(range(d_interp[1], d_interp[end], length(data_x_new)), d2)
itp_data_d1 = itp_d1.(d_interp)
itp_data_d2 = itp_d2.(d_interp)
# analyze peaks
ind_d1, heights_d1 = findmaxima(itp_data_d1)
ind_d1, proms_d1 = peakproms(ind_d1, itp_data_d1)
ind_d1, widths_d1, edges_d1 = peakwidths(ind_d1, itp_data_d1, proms_d1)
widths_d1_a = widths_d1 .* ((d_interp[end]-d_interp[1])/length(itp_data_d1))
for j in eachindex(ind_d1)
    push!(d_peaks_data, (proms_d1[j], widths_d1_a[j]))
end
ind_d2, heights_d2 = findmaxima(itp_data_d2)
ind_d2, proms_d2 = peakproms(ind_d2, itp_data_d2)
ind_d2, widths_d2, edges_d2 = peakwidths(ind_d2, itp_data_d2, proms_d2)
widths_d2_a = widths_d2 .* ((d_interp[end]-d_interp[1])/length(itp_data_d2))
for j in eachindex(ind_d2)
    push!(d_peaks_data, (proms_d2[j], widths_d2_a[j]))
end


# present results
# prominences
av_prom_x = mean([x_peaks_data[j][1] for j in eachindex(x_peaks_data)])
std_prom_x = std([x_peaks_data[j][1] for j in eachindex(x_peaks_data)])

av_prom_y = mean([y_peaks_data[j][1] for j in eachindex(y_peaks_data)])
std_prom_y = std([y_peaks_data[j][1] for j in eachindex(y_peaks_data)])

av_prom_xy = mean([[x_peaks_data[j][1] for j in eachindex(x_peaks_data)]; [y_peaks_data[j][1] for j in eachindex(y_peaks_data)]])
std_prom_xy = std([[x_peaks_data[j][1] for j in eachindex(x_peaks_data)]; [y_peaks_data[j][1] for j in eachindex(y_peaks_data)]])

av_prom_d = mean([d_peaks_data[j][1] for j in eachindex(d_peaks_data)])
std_prom_d = std([d_peaks_data[j][1] for j in eachindex(d_peaks_data)])

# widths in units of unit cell length
av_width_x = mean([x_peaks_data[j][2] for j in eachindex(x_peaks_data)])
std_width_x = std([x_peaks_data[j][2] for j in eachindex(x_peaks_data)])

av_width_y = mean([y_peaks_data[j][2] for j in eachindex(y_peaks_data)])
std_width_y = std([y_peaks_data[j][2] for j in eachindex(y_peaks_data)])

av_width_xy = mean([[x_peaks_data[j][2] for j in eachindex(x_peaks_data)]; [y_peaks_data[j][2] for j in eachindex(y_peaks_data)]])
std_width_xy = std([[x_peaks_data[j][2] for j in eachindex(x_peaks_data)]; [y_peaks_data[j][2] for j in eachindex(y_peaks_data)]])

av_width_d = mean([d_peaks_data[j][2] for j in eachindex(d_peaks_data)])
std_width_d = std([d_peaks_data[j][2] for j in eachindex(d_peaks_data)])