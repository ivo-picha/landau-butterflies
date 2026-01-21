# contains functions related to Landau level basis and wannierization
module States



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
    end

    hf_prev = norm0 * expfac                 # h₀
    hf_curr = sqrt(2f0) * x * hf_prev        # h₁

    @inbounds for k in 1:n-1
        kf = Float32(k)
        hf_next =
            sqrt(2f0 / (kf + 1f0)) * x * hf_curr -
            sqrt(kf / (kf + 1f0)) * hf_prev

        hf_prev, hf_curr = hf_curr, hf_next
    end

    return hf_curr
end


# basis state; Fourier transform of a Landau level in the Landau gauge
function basis_state(x::Float32,y::Float32,n::Int,m::Int,X::Number,Y::Number,p::Int,q::Int,a::Number)
    Kymax = 5 # cutoff in Ky summation; doesn't need to be very large if you are looking at first unit cell

    phi = Float32(p/q)
    lb = Float32(a * sqrt(2π/phi))
    ky(Ky) = Float32(X/(q*a) + 2π*p*Ky/a + 2π*m/a)
    xi(Ky) = x/lb - ky(Ky)*lb

    bs = 0f0
    for Ky in -Kymax:Kymax
        ss = hermite_function(n, xi(Ky))
        ss *= exp(im*(ky(Ky)*y - Ky*Y))
        bs += ss
    end
    bs *= Float32(1/(2π*sqrt(lb*a)))
    # add also normalization over Ky? don't think so
    return bs
end


# get the full eigenstate wavefunction at (x,y) given the eigenvector in LL basis
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


















end