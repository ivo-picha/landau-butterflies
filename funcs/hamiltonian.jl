module Hamil
using LinearAlgebra

# constants
const ħ = 6.62607015e-34/(2π);  # Planck constant [J s]
const e = 1.602176634e-19;      # elementary charge [C]
const m_e = 9.1093837139e-31;   # electron mass [kg];


# recursive laguerre that doesn't need large memory allocation
function laguerre_recursive(α::Number, n::Int64, x::Number)
    if n == 0
        return 1.0
    elseif n == 1
        return -x + α + 1
    else
        Lnm2 = 1.0
        Lnm1 = -x + α + 1
        Ln = 0.0
        for k in 2:n
            Ln = ((2k - 1 + α - x) * Lnm1 - (k - 1 + α) * Lnm2) / k
            Lnm2 = Lnm1
            Lnm1 = Ln
        end
        return Ln
    end
end

# laguerre with the pi/phi ^ n-m/2 factor in the recursion
function laguerre_scaled(n::Int64, m::Int64, x::Number)
    k = min(n,m)
    K = max(n,m)
    d = K - k
    if k == 0
        return x^(d/2)
    elseif k == 1
        return (-x + d + 1) * x^(d/2)
    else
        Ljm2 = x^(d/2)
        Ljm1 = (-x + d + 1) * x^(d/2)
        Lj = 0.0
        for j in 2:k
            Lj = ((2j - 1 + d - x) * Ljm1 - (j - 1 + d) * Ljm2) / j
            Ljm2 = Ljm1
            Ljm1 = Lj
        end
        return Lj
    end
end

function sqrt_factorial_ratio(n::Int64, m::Int64)
    r = 1.0
    for k in (min(n,m)+1):max(n,m)
        r /= sqrt(k)
    end
    return r
end

# LOG SCALE LAGUERRE function
function log_laguerre_scaled(n::Int, m::Int, x::Number)
    k = min(n, m)
    K = max(n, m)
    d = K - k

    logx = log(x)
    log_prefactor = (d / 2) * logx

    if k == 0
        return log_prefactor, 1.0
    end

    if k == 1
        val = -x + d + 1
        return log_prefactor + log(abs(val)), sign(val)
    end

    # Store (log|L|, sign(L))
    logLjm2, signLjm2 = 0.0, 1.0          # corresponds to L₀ = 1
    logLjm1, signLjm1 = log(abs(-x + d + 1)), sign(-x + d + 1)

    for j in 2:k
        a = (2j - 1 + d - x)
        b = (j - 1 + d)

        # Compute term1 = a * Ljm1
        logt1 = log(abs(a)) + logLjm1
        signt1 = sign(a) * signLjm1

        # Compute term2 = b * Ljm2
        logt2 = log(b) + logLjm2
        signt2 = signLjm2

        # term = term1 - term2
        if logt1 > logt2
            logLj = logt1 + log1p(-signt2 * signt1 * exp(logt2 - logt1))
            signLj = signt1
        else
            logLj = logt2 + log1p(-signt1 * signt2 * exp(logt1 - logt2))
            signLj = -signt2
        end

        # divide by j
        logLj -= log(j)

        logLjm2, signLjm2 = logLjm1, signLjm1
        logLjm1, signLjm1 = logLj, signLj
    end

    return log_prefactor + logLjm1, signLjm1
end


# matrix elements (Θ in overleaf)
#T(n::Int64, m::Int64, phi::Float32) = Float32(exp(-Float32(π)/(Float32(2.0)*phi)) * sqrt(bigfac(min(n,m))/bigfac(max(n,m))) * (Float32(π)/phi)^(abs(n-m)/2) * mylaguerre2big(abs(n-m),min(n,m),Float32(π)/phi))
function T(n::Int64, m::Int64, phi::Float32)
    tt = exp(-Float32(π)/(Float32(2.0)*phi))
    tt *= sqrt_factorial_ratio(n,m)
    lls = log_laguerre_scaled(n, m, Float32(π)/phi)
    tt *= exp(lls[1]) * lls[2]
    return Float32(tt)
end

# LL energy in eV
E_LL(n::Int64, phi::Float32, a::Float32) = Float32((n + 0.5) * ħ^2 * 2π * phi / (e * m_e * a^2))

# work in basis (n,m) that goes as (n=0,m=0), (0,1), ..., (0,(p-1)), (1,0),.....
# hamiltonian is composed of A(on diagional) and B matrices that are nonzero on the 3 diagonals

# block matrices on the diagonal -- correspond to (n,k)*(n,k') elements
function matA(n::Int64, phi::Float32, X::Float32, Y::Float32, U0::Float32, a::Float32, p::Int64)
    # get elements from above for the n-th LL
    TnnU = T(n, n, phi) * U0
    En = E_LL(n, phi, a)

    # diagonal elements
    d(j) = En + TnnU * cos(X/Float32(p) + j* Float32(2π)/phi)

    # elements on the two second diagonals
    f = exp(-im*Y/Float32(p)) * TnnU / Float32(2.0)
    fc = conj(f)

    #construct matrix
    mat::Matrix{ComplexF32} = diagm(0 => [d(j) for j = 0:(p-1)], 1 => f * ones(p-1), -1 => fc * ones(p-1))
    mat[1, p] += fc
    mat[p, 1] += f

    return mat
end

# matrix elements away from diagonal; LL mixing
function matB(n::Int64, m::Int64, phi::Float32, X::Float32, Y::Float32, U0::Float32, p::Int64)
    # get elements from above for LL n and m
    TnmU = T(n,m,phi) * U0

    # diagonal elements
    d(j) = TnmU * cos(X/Float32(p) + j*Float32(2π)/phi + abs(m-n)*Float32(π)/Float32(2.0))

    # elements on the two second diagonals
    
    f = exp(-im*Y/Float32(p)) * TnmU / Float32(2.0) * (-1)^((m-n+abs(m-n))/2)     # upper off diagonal
    fc = conj(f) * (-1)^(m-n)                                                     # lower off diagonal

    # f = exp(-im*Y/Float32(p)) * TnmU / Float32(2.0)*(-1)^(m-n)      # upper off diagonal
    # fc = conj(f) * (-1)^(m-n)                                       # lower off diagonal

    #construct matrix
    mat::Matrix{ComplexF32} = diagm(0 => [d(j) for j = 0:(p-1)], 1 => f * ones(p-1), -1 => fc * ones(p-1))
    mat[1, p] += fc
    mat[p, 1] += f

    return mat
end

# row n of the hamiltonian matrix
function get_ham_row(n::Int64, phi::Float32, X::Float32, Y::Float32, U0::Float32, a::Float32, p::Int64, NLL::Int64)
    Z = zeros(ComplexF32, p, p) # zero matrix for filling in empty blocks on lower triangle; autofilled by Hermitian()
    if n == 0
        arrayB = reduce(hcat, [matB(n, m, phi, X, Y, U0, p) for m = Int64(n+1):Int64(NLL)])
        return [matA(n, phi, X, Y, U0, a, p) arrayB]
    elseif n == NLL
        arrayZ = reduce(hcat, [Z for m = Int64(0):Int64(n-1)])
        return [arrayZ matA(n, phi, X, Y, U0, a, p)]
    else
        arrayZ = reduce(hcat, [Z for m = Int64(0):Int64(n-1)])
        arrayB = reduce(hcat, [matB(n, m, phi, X, Y, U0, p) for m = Int64(n+1):Int64(NLL)])
        return [arrayZ matA(n, phi, X, Y, U0, a, p) arrayB]
    end    
end

# full hamiltonian matrix
function get_full_ham(phi::Float32, X::Float32, Y::Float32, U0::Float32, a::Float32, p::Int64, NLL::Int64)
    if NLL == 0 # special case of only lowest LL
        ham = matA(0, phi, X, Y, U0, a, p)
    else
        rows_vector = [get_ham_row(n, phi, X, Y, U0, a, p, NLL) for n = Int64(0):Int64(NLL)]
        ham = reduce(vcat, rows_vector)
    end
    return Hermitian(ham)
end



end