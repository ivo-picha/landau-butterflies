module Wannierize

using Plots
using LinearAlgebra

include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil  


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

    bs = 0f0 + im*0f0
    for Ky in -Kymax:Kymax
        ky = Float32(X/(q*a) + 2π*p*Ky/a + 2π*m/a)
        xi = x/lb - ky*lb
        ss = hermite_function(n, xi)
        ss *= exp(im*(ky*y - Ky*Y)) # minus or plus?! FT determined; minus gives correct result and Chern number
        bs += ss
    end
    bs *= Float32(1/(2π*sqrt(lb*a)))
    # add also normalization over Ky? don't think so
    return bs
end


# get the full eigenstate wavefunction at (x,y) given the eigenvector in LL basis
function get_eigenstate(x::Float32, y::Float32, X::Real, Y::Real, evec::Vector{ComplexF32}, p::Int, q::Int, a::Number, LLmax::Int)::ComplexF32
    
    if p*(LLmax+1) != length(evec)
        println("ERROR: Unexpected eigenvector size")
        exit(1)
    end

    psi = 0f0 + im*0f0
    for idx in 0:(length(evec)-1)
        n = fld(idx, p)
        m = idx%p
        bs = basis_state(x, y, n, m, X, Y, p, q, a)
        psi += evec[idx+1]*bs
    end
    return psi
end

# populate a whole array and normalize over a unit cell (0,qa),(0,a)
function get_bloch_array(xrange, yrange, X::Real, Y::Real, evec::Vector{ComplexF32}, p::Int, q::Int, a::Number, LLmax::Int, dx, dy)
    array = Array{ComplexF32}(undef, length(xrange), length(yrange))
    norm = 0f0
    for (ix,x) in enumerate(xrange)
        for (jy,y) in enumerate(yrange)
            wf_xy = get_eigenstate(x,y,X,Y,evec,p,q,a,LLmax)
            array[ix,jy] = wf_xy

            if 0 <= x < q*a && 0<= y < a
                norm += abs2(wf_xy)
            end
        end
    end
    array[:,:] ./= sqrt(norm*dx*dy)
    return array
end



# populate a whole array, peel magnetic phase and normalize over a unit cell (0,qa),(0,a)
function get_bloch_array_magnetic_peel(xrange, yrange, X::Real, Y::Real, evec::Vector{ComplexF32}, p::Int, q::Int, a::Number, LLmax::Int, dx, dy)
    array = Array{ComplexF32}(undef, length(xrange), length(yrange))
    norm = 0f0
    pi32=Float32(π)
    for (ix,x) in enumerate(xrange)
        for (jy,y) in enumerate(yrange)
            wf_xy = get_eigenstate(x,y,X,Y,evec,p,q,a,LLmax)
            array[ix,jy] = wf_xy * exp(-im*2*pi32*p*(y/a)*floor(x/(q*a)))

            if 0 <= x < q*a && 0<= y < a
                norm += abs2(wf_xy)
            end
        end
    end
    array[:,:] ./= sqrt(norm*dx*dy)
    return array
end


# populate a whole array and normalize over a unit cell (0,qa),(0,a)
function get_bloch_periodic(xrange, yrange, X::Real, Y::Real, evec::Vector{ComplexF32}, p::Int, q::Int, a::Number, LLmax::Int,dx, dy)
    array = Array{ComplexF32}(undef, length(xrange), length(yrange))
    norm = 0f0
    kx=Y/(q*a)
    ky=X/(q*a)
    pi32=Float32(π)
    for (ix,x) in enumerate(xrange)
        for (jy,y) in enumerate(yrange)
            wf_xy = get_eigenstate(x,y,X,Y,evec,p,q,a,LLmax)
            array[ix,jy] = wf_xy * exp(-im*2*pi32*p*(y/a)*floor(x/(q*a))) * exp(-im*ky*y) * exp(-im*kx*x)

            if 0 <= x < q*a && 0<= y < a
                norm += abs2(wf_xy)
            end
        end
    end
    array[:,:] ./= sqrt(norm*dx*dy)
    return array
end


function get_trial_from_u(n::Integer, xrange_muc, yrange_muc,X::Real,Y::Real, evec::Vector{ComplexF32}, p::Int, q::Int, a::Number, LLmax::Int,dx, dy)
    trial = zeros(ComplexF32, length(xrange_muc), length(yrange_muc))
    steps_uc = Int(length(xrange_muc)/q)
    xrange_uc = range((n-1)*a, n*a, length=steps_uc)
    trial[(n-1)*steps_uc+1:n*steps_uc,:] = get_bloch_array(xrange_uc,yrange_muc,X,Y,evec,p,q,a,LLmax,dx, dy)

    #normalize
    norm = sqrt(sum(abs2.(trial)) * dx*dy)
    trial[:,:] ./= norm

    return trial
end
    

# plot an array 
function plot_array(xrange,yrange,array,a,q)
    array[:,:] ./= maximum(abs.(array[:,:])) # normalize for plotting
    plt_abs = heatmap(xrange./a, yrange./a, abs.(array[:,:])', color=:viridis, aspect_ratio=1)
    plt_arg = heatmap(xrange./a, yrange./a, angle.(array[:,:])', color=:hsv, aspect_ratio=1)
    pltcombi = plot(plt_abs, plt_arg, layout = (2,1), size = (q*500,1000), framestyle=:box)
    return pltcombi
end


function get_U_matrix_and_energies(p,q,X,Y,U0,a,LLmax,Nuc,xrange_muc,yrange_muc,dx,dy,trialwf)

    phi = Float32(p/q)

    H = Hamil.get_full_ham(phi, Float32(X), Float32(Y), U0, a, p, LLmax)
    evals, evecs = eigen(H)


    bloch_wf_array_muc = Array{ComplexF32}(undef, length(xrange_muc), length(yrange_muc), q)
    for m=1:q
        bloch_wf_array_muc[:,:,m] = Wannierize.get_bloch_array(xrange_muc,yrange_muc,X,Y,evecs[:,m],p,q,a,LLmax,dx,dy)
    end

    Amat = zeros(ComplexF32, q, q)

    for n=1:q
        for m=1:q
            Amat[m,n] = sum( conj.(bloch_wf_array_muc[:,:,m]) .* trialwf[:,:,n] ) * dx*dy
        end
    end

    Smat = Hermitian(Amat' * Amat)
    Umat = Amat * sqrt(inv(Smat))

    if norm(Umat' * Umat - I) > 0.01
        println("WARNING: U matrix is not unitary!")
    end

    return Umat, evals
    
end


function get_loewdin_array(p,q,X,Y,U0,a,LLmax,Nuc,xrange,yrange,xrange_muc,yrange_muc,dx,dy,trialwf)

    phi = Float32(p/q)

    H = Hamil.get_full_ham(phi, Float32(X), Float32(Y), U0, a, p, LLmax)
    evals, evecs = eigen(H)

    Umat, ens = get_U_matrix_and_energies(p,q,X,Y,U0,a,LLmax,Nuc,xrange_muc,yrange_muc,dx,dy,trialwf)

    bloch_wf_array = Array{ComplexF32}(undef, length(xrange), length(yrange), q)
    for m=1:q
        bloch_wf_array[:,:,m] = Wannierize.get_bloch_array(xrange,yrange,X,Y,evecs[:,m],p,q,a,LLmax,dx,dy)
    end


    loewdin_wf_array = zeros(ComplexF32, length(xrange), length(yrange), q)
    for n=1:q
        for m=1:q
            loewdin_wf_array[:,:,n] .+= Umat[m,n] .* bloch_wf_array[:,:,m]
        end
    end

    return loewdin_wf_array
    
end


# Landau gauge ground state of harmonic + magnetic; adds phase Bxy
function gaussian_LG(x::Float32,y::Float32,Rx::Float32,Ry::Float32, phi::Float32, U0::Float32, a::Float32)
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
    return exp(-((x-Rx)^2 + (y-Ry)^2)/(2*σsq)) #* exp(im*Float32(2π)*phi*x*y/(a^2))
end


function write_complex_matrix_in_io(io, M::AbstractMatrix{<:Complex}; digits::Int=7)
    nrows, ncols = size(M)

    for i in 1:nrows
        print(io, "  ")

        for j in 1:ncols
            z = M[i, j]

            r = abs(z)
            θπ = angle(z) / π

            r_str = string(round(r; digits=digits))
            θπ_str = string(round(θπ; digits=digits))

            entry = "$(r_str) ∠ $(θπ_str)π"

            print(io, lpad(entry, 22))

            if j < ncols
                print(io, "   ")
            end
        end

        if i < nrows
            println(io, ";")
        else
            println(io)
        end
    end

end


function get_bloch_Rxmax(xrange,array,a)
    xyc = findmax(abs2, array)
    x = xrange[xyc[2][1]]/a
    return round(x; digits=2)
end


function get_all_trial_wavefunctions(p,q,U0,a,LLmax,xrange_muc,yrange_muc,dx,dy)
    phi = Float32(p/q)
    H0 = Hamil.get_full_ham(phi, 0f0, 0f0, U0, a, p, LLmax)
    evals0, evecs0 = eigen(H0)
    trialwfs = zeros(ComplexF32, length(xrange_muc), length(yrange_muc),q)

    steps_uc = Int(length(xrange_muc)/q)

    for m=1:q
        bloch_array = get_bloch_array(xrange_muc,yrange_muc,0f0,0f0,evecs0[:,m],p,q,a,LLmax,dx,dy)
        Rxmax = get_bloch_Rxmax(xrange_muc,bloch_array,a)

        xrange_uc = range(Int(round(Rxmax-0.5))*a, Int(round(Rxmax+0.5))*a, length=steps_uc)

        trialwfs[Int(round(Rxmax-0.5))*steps_uc+1:Int(round(Rxmax+0.5))*steps_uc,:,m] = get_bloch_array(xrange_uc,yrange_muc,0f0,0f0,evecs0[:,m],p,q,a,LLmax,dx, dy)
        norm = sqrt(sum(abs2.(trialwfs[Int(round(Rxmax-0.5))*steps_uc+1:Int(round(Rxmax+0.5))*steps_uc,:,m])) * dx*dy)
        trialwfs[:,:,m] ./= norm
    end

    return trialwfs

end


end