# contains functions related to Landau level basis and wannierization
module States

using Plots, Measures

function fermi_dirac(en::Number, eF::Number, eT::Number)
    return 1/(exp((en-eF)/eT)+1)
end

# recursive hermite polynomial; normalized
function hermite_function(n::Int, x::Float32)::Float32
    # common factor
    expfac = exp(-x*x/2f0)
    norm0  = π^(-1f0/4f0)

    if n == 0
        return norm0 * expfac
    elseif n == 1
        return sqrt(2f0) * x * norm0 * expfac
    else

        hf_prev = norm0 * expfac                 # h₀
        hf_curr = sqrt(2f0) * x * hf_prev        # h₁

        for k in 1:n-1
            kf = Float32(k)
            hf_next =
                sqrt(2f0 / (kf + 1f0)) * x * hf_curr -
                sqrt(kf / (kf + 1f0)) * hf_prev

            hf_prev, hf_curr = hf_curr, hf_next
        end

        return hf_curr
    end
end


# basis state; Fourier transform of a Landau level in the Landau gauge
function basis_state(x::Float32,y::Float32,n::Int,m::Int,X::Number,Y::Number,p::Int,q::Int,a::Number)
    Kymax = 5 # cutoff in Ky summation; doesn't need to be very large if you are looking at first unit cell

    phi = Float32(p/q)
    lb = Float32(a / sqrt(2π*phi))

    bs = 0f0
    for Ky in -Kymax:Kymax
        ky = Float32(X/(q*a) + 2π*p*Ky/a + 2π*m/a)
        xi = x/lb - ky*lb
        ss = hermite_function(n, xi)
        ss *= exp(im*(ky*y - Ky*Y))
        bs += ss
    end
    bs *= Float32(1/(2π*sqrt(lb*a)))
    # add also normalization over Ky? don't think so
    return bs
end


# get the full eigenstate wavefunction at (x,y) given the eigenvector in LL basis; for phi 1 use
function get_eigenstate(x::Float32, y::Float32, state::Tuple{Float32, Vector{ComplexF32}, Float32, Float32}, p::Int, q::Int, a::Number, LLmax::Int)::ComplexF32
    en, evec, X, Y = state
    n_states = length(evec)
    n_LL = LLmax + 1
    if n_LL*p != n_states
        error("Number of states in eigenvector does not match number of Landau levels times p.")
    end
    psi = 0f0 + im*0f0
    for idx in 0:n_states-1
        n = fld(idx, p)
        m = idx%p
        bs = basis_state(x, y, n, m, X, Y, p, q, a)
        psi += evec[idx+1]*bs
    end
    return psi
end

# get the full eigenstate wavefunction at (x,y) given the eigenvector in LL basis; for general use
function get_eigenstate_XY(x::Float32, y::Float32, state::Tuple{Float32, Vector{ComplexF32}}, X::Float32, Y::Float32, p::Int, q::Int, a::Number, LLmax::Int)::ComplexF32
    en, evec = state
    n_states = length(evec)
    n_LL = LLmax + 1
    if n_LL*p != n_states
        error("Number of states in eigenvector does not match number of Landau levels times p.")
    end
    psi = 0f0 + im*0f0
    for idx in 0:n_states-1
        n = fld(idx, p)
        m = idx%p
        bs = basis_state(x, y, n, m, X, Y, p, q, a)
        psi += evec[idx+1]*bs
    end
    return psi
end


# gaussian trial orbital wavefunction for wannierization
function gaussian(x::Float32,y::Float32,X::Float32,Y::Float32,σ::Float32)
    # gaussian trial orbital wavefunction for wannierization
    return exp(-((x-X)^2 + (y-Y)^2)/(2*σ^2))
end

# Landau gauge ground state of harmonic + magnetic; adds phase Bxy
function gaussian_LG(x::Float32,y::Float32,X::Float32,Y::Float32, phi::Float32, U0::Float32, a::Float32)
    phi = Float64(phi)
    U0 = Float64(U0)
    a = Float64(a)
    h = 6.62607015e-34;  # Planck constant [J s]
    e = 1.602176634e-19;      # elementary charge [C]
    m_e = 9.1093837139e-31;   # electron mass [kg]
    # gaussian trial orbital wavefunction for wannierization
    ωhsq = 4.0 * π^2 * abs(U0) * e / (m_e * a^2) # harmonic frequency squared
    ωcsq4 = phi^2 * h^2 / (16.0 * m_e^2 * a^4)
    σsq = h / (2π * m_e * sqrt(ωhsq + ωcsq4))
    return exp(-((x-X)^2 + (y-Y)^2)/(2*σsq)) * exp(-im*Float32(π)*phi*x*y/(a^2))
end


# calculate wannier function at (x,y) and centre R
function get_wannier_array(q::Int, Rx::Float32, Ry::Float32, loewdin_array::Array{ComplexF32,4}, X_list::Vector{Float32}, Y_list::Vector{Float32})
    wannier_array = Array{ComplexF32}(undef, size(loewdin_array,1), size(loewdin_array,2))
    for i in 1:size(loewdin_array,1)
        for j in 1:size(loewdin_array,2)
            ws = 0f0 + im*0f0
            for (Xidx,X) in enumerate(X_list)
                for (Yidx,Y) in enumerate(Y_list)
                    ws += loewdin_array[i,j,Xidx,Yidx] * exp(im*(-X*Ry + Y*Rx)/q)
                end
            end
            wannier_array[i,j] = ws / (length(X_list)*length(Y_list))
        end
    end
    norm = sum(abs.(wannier_array))
    return wannier_array ./ norm # still needs normalization a^2/(2pi)^2 or similar
end


# plot a specific wannier function, given its array representation
function plot_wannier(wannier_array::Array{ComplexF32,2}, xplotrange::Vector{Float32}, yplotrange::Vector{Float32}, out_folder::String, U0::Float32, a::Float32, p::Int, q::Int, Rx::Float32, Ry::Float32, Ngrid::Int)
    # visualize wannier function
    wf_real = real.(wannier_array)
    wf_imag = imag.(wannier_array)
    wf_abs2 = abs2.(wannier_array)
    wf_arg = angle.(wannier_array)
    plt_r = heatmap(xplotrange./a, yplotrange./a, wf_real', xlabel="x/a", ylabel="y/a", title="Re(wᵣ)", aspect_ratio=1, size=(500,500), legend = false);
    plt_i = heatmap(xplotrange./a, yplotrange./a, wf_imag', xlabel="x/a", ylabel="y/a", title="Im(wᵣ)", aspect_ratio=1, size=(500,500), legend = false);
    plt_a = heatmap(xplotrange./a, yplotrange./a, sqrt.(wf_abs2)', xlabel="x/a", ylabel="y/a", title="|wᵣ|", aspect_ratio=1, size=(500,500), legend = false);
    plt_th = heatmap(xplotrange./a, yplotrange./a, reshape(wf_arg', (Ngrid, Ngrid))', xlabel="x/a", ylabel="y/a", title="arg(wᵣ)", aspect_ratio=1, size=(500,500), legend = false);
    plt_combo = plot(plt_r, plt_i, plt_a, plt_th, layout = grid(2, 2, hgap = 2mm, vgap = 2mm),
        aspect_ratio = 1,
        size = (800,800),
        margin = 1mm, plot_title="wᵣ @ rx/a=$Rx, ry/a=$Ry; ϕ=$(p)/$(q), U0=$(U0) eV"
    )
    savefig(plt_combo, joinpath(out_folder, "wannier_function_$(U0)_phi$(p)_$(q)_rx$Rx-ry$Ry.png"))
    println("Saved wannier function plot to $(joinpath(out_folder, "wannier_function_$(U0)_phi$(p)_$(q)_rx$Rx-ry$Ry.png"))")
end


function get_derivative_x(array::Array{ComplexF64,2}, dx::Number)::Array{ComplexF64,2}
    Ngridx, Ngridy = size(array)
    deriv_array = Array{ComplexF64}(undef, Ngridx, Ngridy)
    for j in 1:Ngridy
        for i in 1:Ngridx
            if i == 1
                deriv_array[i,j] = (array[i+1,j] - array[i,j]) / dx
            elseif i == Ngridx
                deriv_array[i,j] = (array[i,j] - array[i-1,j]) / dx
            else
                deriv_array[i,j] = (array[i+1,j] - array[i-1,j]) / (2f0*dx)
            end
        end
    end
    return deriv_array
end

function get_derivative_y(array::Array{ComplexF64,2}, dy::Number)::Array{ComplexF64,2}
    Ngridx, Ngridy = size(array)
    deriv_array = Array{ComplexF64}(undef, Ngridx, Ngridy)
    for i in 1:Ngridx
        for j in 1:Ngridy
            if j == 1
                deriv_array[i,j] = (array[i,j+1] - array[i,j]) / dy
            elseif j == Ngridy
                deriv_array[i,j] = (array[i,j] - array[i,j-1]) / dy
            else
                deriv_array[i,j] = (array[i,j+1] - array[i,j-1]) / (2f0*dy)
            end
        end
    end
    return deriv_array
end

# build hamiltonian, acting on wannier state:
# T w = hbar^2/2m *(-dx^2 - dy^2 + 2ei/hbar * Bx dy + e^2/hbar^2 B^2 x^2 ) w
# eB/ħ = ϕ /(2π a^2)
function H_on_wannier(wannier_array::Array{ComplexF32,2}, xrange::Vector{Float32}, yrange::Vector{Float32}, phi::Float32, U0::Float32, a::Float32)::Array{ComplexF32}
    
    # work in 64 bit precision
    phi = Float64(phi)
    U0 = Float64(U0)
    a = Float64(a)
    xrange = Float64.(xrange)
    yrange = Float64.(yrange)
    wannier_array = ComplexF64.(wannier_array)

    ħ = 6.62607015e-34/(2π);  # Planck constant [J s]
    e = 1.602176634e-19;      # elementary charge [C]
    m_e = 9.1093837139e-31;   # electron mass [kg]

    dx = xrange[2] - xrange[1]
    dy = yrange[2] - yrange[1]
    Ngridx, Ngridy = size(wannier_array)

    # init empty array; add term by term, then multiply
    h_wannier = zeros(ComplexF64, Ngridx, Ngridy)

    #derivatives
    wdx = get_derivative_x(wannier_array, dx)
    wdx2 = get_derivative_x(wdx, dx)
    wdy = get_derivative_y(wannier_array, dy)
    wdy2 = get_derivative_y(wdy, dy)

    # arrays of x values
    x_array = Array{Float64}(undef, Ngridx, Ngridy)
    for i = 1:Ngridx
        for j = 1:Ngridy
            x_array[i,j] = xrange[i]
        end
    end

    y_array = Array{Float64}(undef, Ngridx, Ngridy)
    for i = 1:Ngridx
        for j = 1:Ngridy
            y_array[i,j] = yrange[j]
        end
    end

    # array of x^2 values
    x2_array = x_array.^2

    # eB/ħ
    ebh = Float32(2π * phi / (a^2)) # 2pi up or down?

    h_wannier += -1.0.*wdx2
    h_wannier += -1.0.*wdy2
    h_wannier += 2.0*im*ebh.* (x_array .* wdy)
    h_wannier += ebh^2 .* (x2_array .* wannier_array)
    h_wannier *= ħ^2 /(2 * m_e)
    h_wannier ./= e # divide by e to get energy in eV

    # add potential term
    pot = cos.(x_array.*Float64(2π/a)) 
    pot += cos.(y_array.*Float64(2π/a))
    pot = ComplexF64.(pot)
    pot .*= wannier_array
    pot *= U0

    h_wannier += pot

    return h_wannier 
end




end