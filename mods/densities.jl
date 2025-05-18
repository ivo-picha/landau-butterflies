module Dens

using LinearAlgebra
using Dierckx
using Peaks
using Interpolations
using Statistics: mean, std, cor

# Fermi-Dirac distribution
function fermi_dirac(en::Float64, eF::Float64, eT::Float64)
    return 1/(exp((en-eF)/eT)+1)
end

# cut off with a bonus so that spectrum can be smeared later without a sharp cutoff
function discard_high_energies_smear(states_vec_sorted::Vector{Tuple{Float64, Float64, Float64, Vector{ComplexF64}}}, phi::Float64, NLL::Int64, np::Float64, TeV::Float64)
    N_en = length(states_vec_sorted)
    np_max = phi * (NLL+1)
    N_cut = round(Int, N_en * np / np_max)
    # find position on which energy is suff large
    if N_cut >= N_en
        println("$np particles per unit cell is more than the allowed maximum ($np_max) at given flux, NLL. Taking $np_max particles per unit cell instead. Add more LLs to the calculation!")
        return (states_vec_sorted, states_vec_sorted[end][1])
    else
        EF = states_vec_sorted[N_cut][1] # Fermi energy
        E_cutoff = EF + 4*TeV
        energies_sorted = first.(states_vec_sorted)
        N_cut_plus = searchsortedlast(energies_sorted, E_cutoff; lt = <)
        if N_cut_plus >= N_en
            println("Temperature smear clips the end of the spectrum. Add more LLs to the calculation!")
            return (states_vec_sorted, EF)
        else
            return (states_vec_sorted[1:N_cut_plus], EF) # returns the list of states + fermi energy 
        end
    end
end


# recursively define a hermite polynomial
function hermite_r(x::Float64, n::Int)::Float64
    if n == 0
        return 1.0
    elseif n == 1
        return 2*x
    else
        return 2*x* hermite_r(x, n-1) - 2*(n-1)* hermite_r(x, n-2)
    end
end

# orthogonal hermite function
function hermite_function(x::Float64, n::Int64)
    An = 1/sqrt(2^n * factorial(n) * sqrt(π))
    return An * hermite_r(x,n) * exp(-x^2 /2)
end

function landau_lvl_wf(x::Float64, y::Float64, n::Int, m::Int, ky0::Float64, Ky::Int, phi::Float64, a::Float64, p::Int)::ComplexF64
    lB = a / sqrt(2π * phi)
    ky = ky0 + 2π*m/a + 2π*p*Ky/a
    xix = x/lB - ky*lB
    return 1/a * hermite_function(xix, n) * exp(im * ky * y)
end


function get_density_list(xyplotlist::Vector{Tuple{Float64,Float64}}, state::Tuple{Float64, Float64, Float64, Vector{ComplexF64}}, nm_list::Vector{Tuple{Int,Int}}, phi::Float64, a::Float64, p::Int, NKy::Int=2)
    dens_list = Float64[];
    for xy in xyplotlist
        x = xy[1]
        y = xy[2]
        wf = sum([state[4][j] * exp(-im*Ky*state[3]) * landau_lvl_wf(x,y,nm_list[j][1],nm_list[j][2],state[2],Int(Ky),phi,a,p) for j in eachindex(nm_list) for Ky = -NKy:NKy])
        push!(dens_list, real(wf*conj(wf)))
    end
    return dens_list
end


# interpolate the electronic density data
function square_interpolate_d(data_x::Vector{Float64}, data_y::Vector{Float64}, data_z::Matrix{Float64}, Nnew::Int=256)
    
    spline_d = Spline2D(data_x, data_y, copy(transpose(data_z)))

    data_x_new = range(data_x[1], data_x[end], Nnew)
    data_y_new = range(data_y[1], data_y[end], Nnew)

    data_z_new = evalgrid(spline_d, data_x_new, data_y_new)

    return [collect(data_x_new), collect(data_y_new), data_z_new]
end




# FUNCTIONS THAT ANALYZE THE ELECTRONIC DENSITY

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


# USING PEAKS.JL
# exract the second main diagonal of a matrix
function antidiag(A::Matrix)
    rows, cols = size(A)
    len = min(rows, cols)
    return [A[i, cols - i + 1] for i in 1:len]
end


# analyze the density data and output prominence and width mean and standard deviation in xy or diagonal
function get_proms_widths(data_x::Vector{Float64}, data_y::Vector{Float64}, data_z::Matrix{Float64}, Nuc::Int)
    fm = findmax(data_z)
    translation_index = round(Int, length(data_x)/Nuc)

    # set the number of steps for an interpolation of the cross-sections
    Nitp = 512
    
    # vertical cross-sections
    x_peaks_data = [];
    x_interp = range(data_x[1], data_x[end], Nitp)
    for i = 1:Nuc
        cs = data_z[Int(mod( fm[2][1] + translation_index*i, translation_index*Nuc)),:]
        itp_i = Interpolations.CubicSplineInterpolation(range(data_x[1], data_x[end], length(data_x)), cs)
        itp_data = itp_i.(x_interp)
    
        #Peaks.jl for analysing cross-sections
        ind_i, heights_i = findmaxima(itp_data)
        ind_i, proms_i = peakproms(ind_i, itp_data)
        ind_i, widths_i, edges_i = peakwidths(ind_i, itp_data, proms_i)
    
        widths_i_a = widths_i .* ((data_x[end]-data_x[1])/length(itp_data))
        for j in eachindex(ind_i)
            push!(x_peaks_data, (proms_i[j], widths_i_a[j]))
        end
    end
    
    # vertical cross-sections
    y_peaks_data = [];
    y_interp = range(data_y[1], data_y[end], Nitp)
    for i = 1:Nuc
        cs = data_z[:, Int(mod( fm[2][2] + translation_index*i, translation_index*Nuc))]
        itp_i = Interpolations.CubicSplineInterpolation(range(data_y[1], data_y[end], length(data_y)), cs)
        itp_data = itp_i.(y_interp)
    
        #Peaks.jl for analysing cross-sections
        ind_i, heights_i = findmaxima(itp_data)
        ind_i, proms_i = peakproms(ind_i, itp_data)
        ind_i, widths_i, edges_i = peakwidths(ind_i, itp_data, proms_i)
    
        widths_i_a = widths_i .* ((data_y[end]-data_y[1])/length(itp_data))
        for j in eachindex(ind_i)
            push!(y_peaks_data, (proms_i[j], widths_i_a[j]))
        end
    end
    
    
    # diagonal cross-sections on the two main diagonals 
    d_peaks_data = [];
    d_interp = range(sqrt(data_x[1]^2 + data_y[1]^2), sqrt(data_x[end]^2 + data_y[end]^2), Nitp)
    d1 = diag(data_z)
    d2 = antidiag(data_z)
    itp_d1 = Interpolations.CubicSplineInterpolation(range(d_interp[1], d_interp[end], length(data_x)), d1)
    itp_d2 = Interpolations.CubicSplineInterpolation(range(d_interp[1], d_interp[end], length(data_x)), d2)
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
    # av_prom_x = mean([x_peaks_data[j][1] for j in eachindex(x_peaks_data)])
    # std_prom_x = std([x_peaks_data[j][1] for j in eachindex(x_peaks_data)])
    
    # av_prom_y = mean([y_peaks_data[j][1] for j in eachindex(y_peaks_data)])
    # std_prom_y = std([y_peaks_data[j][1] for j in eachindex(y_peaks_data)])
    
    av_prom_xy = mean([[x_peaks_data[j][1] for j in eachindex(x_peaks_data)]; [y_peaks_data[j][1] for j in eachindex(y_peaks_data)]])
    std_prom_xy = std([[x_peaks_data[j][1] for j in eachindex(x_peaks_data)]; [y_peaks_data[j][1] for j in eachindex(y_peaks_data)]])
    
    av_prom_d = mean([d_peaks_data[j][1] for j in eachindex(d_peaks_data)])
    std_prom_d = std([d_peaks_data[j][1] for j in eachindex(d_peaks_data)])
    
    # widths in units of unit cell length
    # av_width_x = mean([x_peaks_data[j][2] for j in eachindex(x_peaks_data)])
    # std_width_x = std([x_peaks_data[j][2] for j in eachindex(x_peaks_data)])
    
    # av_width_y = mean([y_peaks_data[j][2] for j in eachindex(y_peaks_data)])
    # std_width_y = std([y_peaks_data[j][2] for j in eachindex(y_peaks_data)])
    
    av_width_xy = mean([[x_peaks_data[j][2] for j in eachindex(x_peaks_data)]; [y_peaks_data[j][2] for j in eachindex(y_peaks_data)]])
    std_width_xy = std([[x_peaks_data[j][2] for j in eachindex(x_peaks_data)]; [y_peaks_data[j][2] for j in eachindex(y_peaks_data)]])
    
    av_width_d = mean([d_peaks_data[j][2] for j in eachindex(d_peaks_data)])
    std_width_d = std([d_peaks_data[j][2] for j in eachindex(d_peaks_data)])

    return (av_prom_xy, std_prom_xy, av_width_xy, std_width_xy, av_prom_d, std_prom_d, av_width_d, std_width_d)
end



end