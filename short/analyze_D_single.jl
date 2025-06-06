
using NPZ
using Dierckx
using Plots
using LinearAlgebra
using Statistics: cor, mean
npz_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/densities_zerofield/ZF_U00.01-a5.0e-9-nu1.2.npz"
plot_save_folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/npzdens"

# interpolate the electronic density data
function square_interpolate_d(data_x::Vector{Float64}, data_y::Vector{Float64}, data_z::Matrix{Float64}, Nnew::Int=256)
    
    spline_d = Spline2D(data_x, data_y, copy(transpose(data_z)))

    data_x_new = range(data_x[1], data_x[end], Nnew)
    data_y_new = range(data_y[1], data_y[end], Nnew)

    data_z_new = evalgrid(spline_d, data_x_new, data_y_new)

    return [collect(data_x_new), collect(data_y_new), data_z_new]
end

a = 50e-10
data = npzread(npz_path)
data_x = data["x"] ./a
data_y = data["y"] ./a
data_z = data["z"]

dxn, dyn, dzn = square_interpolate_d(data_x, data_y, data_z)

heatmap(dxn,dyn,dzn)
# heatmap(data_x,data_y,data_z)

########################## standard dev in R
function get_exp_R(dxn, dyn, dzn)
    norm = sum(dzn)
    Rx = sum([dzn[i,j] * dxn[i] for i in eachindex(dxn) for j in eachindex(dyn)])
    Ry = sum([dzn[i,j] * dyn[j] for i in eachindex(dxn) for j in eachindex(dyn)])
    return (Rx/norm, Ry/norm)
end

function get_exp_Rsq(dxn, dyn, dzn)
    norm = sum(dzn)
    Rsq = sum([dzn[i,j] * (dxn[i]^2 + dyn[j]^2) for i in eachindex(dxn) for j in eachindex(dyn)])
    return Rsq/norm    
end

function get_sigmaR(dxn,dyn,dzn)
    exp_R = get_exp_R(dxn,dyn,dzn)
    return sqrt(get_exp_Rsq(dxn,dyn,dzn)-(exp_R[1]^2 + exp_R[2]^2))
end

########################## prominence
function get_dens_ratio(dzn)
    d_p0, maxind = findmax(dzn)
    d_px1 = dzn[1, maxind[2]]
    d_px2 = dzn[end, maxind[2]]
    d_py1 = dzn[maxind[1], 1]
    d_py2 = dzn[maxind[1], end]
    return 4*d_p0/(d_px1+d_px2+d_py1+d_py2)
end

########################## diamondness
function rotational_symmetry_score(data::Matrix{Float64}; angles=0:5:345)
    rows, cols = size(data)

    (maxval, maxidx) = findmax(data)
    cy, cx = Tuple(maxidx)

    # Convert matrix to centered coordinate grid
    function rotate_coords(y, x, angle_rad, cy, cx)
        # Shift to origin, rotate, shift back
        dy, dx = y - cy, x - cx
        y_rot = cos(angle_rad)*dy - sin(angle_rad)*dx + cy
        x_rot = sin(angle_rad)*dy + cos(angle_rad)*dx + cx
        return round(Int, y_rot), round(Int, x_rot)  # Nearest-neighbor
    end

    # Compute average correlation over rotations
    corrs = Float64[]

    for θ_deg in angles
        θ_rad = θ_deg * π / 180
        rotated = fill(NaN, rows, cols)

        for y in 1:rows, x in 1:cols
            y_r, x_r = rotate_coords(y, x, θ_rad, cy, cx)
            if 1 ≤ y_r ≤ rows && 1 ≤ x_r ≤ cols
                rotated[y, x] = data[y_r, x_r]
            end
        end

        # Compare only valid (non-NaN) entries
        valid_mask = .!isnan.(rotated)
        orig_vals = data[valid_mask]
        rot_vals = rotated[valid_mask]

        push!(corrs, cor(orig_vals, rot_vals))
    end

    return mean(corrs)
end

get_sigmaR(dxn,dyn,dzn)
