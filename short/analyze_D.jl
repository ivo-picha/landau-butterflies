# plot a density plot of the electronic density from a pre-saved .npz file

include(joinpath(dirname(@__DIR__), "mods/params.jl"))
using .Params

using NPZ
using Plots
using Peaks
using Interpolations
using Dierckx

npz_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/local/dens_grids_p1-q1-U0.01-a50.0-N5-n1.0-T10.0.npz"
plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/npzdens"

data = npzread(npz_path)
data_x = data["x"]
data_y = data["y"]
data_z = data["z"]

# data_z_t = data_z
# for i in eachindex(data_x)
#     for j in eachindex(data_y)
#         if i > 16 && j>16
#             data_z_t[i,j] = 0
#         end
#     end
# end

# heatmap(data_x, data_x, data_z_t)

spline = Spline2D(data_x, data_y, copy(transpose(data_z_t)))


Nnew = 128
xgrid_new = range(data_x[1], data_x[end], Nnew)
ygrid_new = range(data_y[1], data_y[end], Nnew)

zgrid_new = evalgrid(spline, xgrid_new, ygrid_new)

heatmap(xgrid_new, ygrid_new, zgrid_new)

Nuc = 2
#function get_d_params(data_x::Vector{Float64}, data_y::Vector{Float64}, data_z::Matrix{Float64}, Nuc = 2)

    fm = findmax(data_z)
    mp = round(Int, length(data_x)/Nuc)

    # perhaps extend to cover all peaks automatically

    data_line_x1 = data_z[fm[2][1],:]
    data_line_x2 = data_z[mod(fm[2][1] + mp, mp*Nuc),:]

    data_line_y1 = data_z[:,fm[2][2]]
    data_line_y2 = data_z[:,mod(fm[2][2]+mp, mp*Nuc)]

    x_interp = range(data_x[1], data_x[end], length(data_x)*10)
    y_interp = range(data_y[1], data_y[end], length(data_y)*10)

    itp_x1 = CubicSplineInterpolation(range(data_x[1], data_x[end], length(data_x)), data_line_x1)
    itp_x2 = CubicSplineInterpolation(range(data_x[1], data_x[end], length(data_x)), data_line_x2)
    itp_y1 = CubicSplineInterpolation(range(data_y[1], data_y[end], length(data_y)), data_line_y1)
    itp_y2 = CubicSplineInterpolation(range(data_y[1], data_y[end], length(data_y)), data_line_y2)

    itpdata_x1 = itp_x1.(x_interp)
    itpdata_x2 = itp_x2.(x_interp)
    itpdata_y1 = itp_y1.(y_interp)
    itpdata_y2 = itp_y2.(y_interp)


    ind_x1, hei_x1 = findmaxima(itpdata_x1)
    ind_x1, proms_x1 = peakproms(ind_x1, itpdata_x1)
    ind_x1, widths_x1, edges_x1... = peakwidths(ind_x1, itpdata_x1, proms_x1)

    widths_x1_a = widths_x1 .* (data_x[end]-data_x[1]/length(itpdata_x1))
    proms_x1

    #plot(x_interp, itpdata_x1)
 
#end

