module Hamil
using Polynomials, SpecialPolynomials
using LinearAlgebra

# constants
const ħ = 6.62607015e-34/(2π);  # Planck constant [J s]
const e = 1.602176634e-19;      # elementary charge [C]
const m_e = 9.1093837139e-31;   # electron mass [kg];

# define an easier to use laguerre polynomial 
function mylaguerre(α::Number, n::Int64, x::Number)
    list1 = push!(zeros(Integer,n),1)
    lag = Laguerre{α}(list1)
    return lag(x)    
end

# define functions which can handle big factorials and powers beyond Int128
function bigfac(n::Int)
    if n > 20
        return factorial(big(n))
    else
        return factorial(n)
    end    
end

function bigbinomial(n::Int, m::Int)
        if n > 30
        return binomial(big(n),big(m))
    else
        return binomial(n,m)
    end    
end

function bigpw2(n::Int)
        if n > 30
        return big(2)^n
    else
        return 2^n
    end    
end

# define an easier to use laguerre polynomial which can handle big numbers without overflow
function mylaguerre2big(α::Number, n::Int64, x::Number)
    return sum([(-1)^k * bigbinomial(n+α,n-k) * big(x)^k /(bigfac(k)) for k=0:n])    
end

# LL energy in eV
E_LL(n::Int64, ξ0::Float64, a::Float64) = (n + 0.5) * ħ^2 / (e * m_e * (ξ0 * a / (2π))^2)

# matrix elements (Θ in overleaf)
Tx(n::Int64, m::Int64, ξ0::Float64) = exp(- ξ0^2 / 4) * (sqrt((bigpw2(n+m)))/(sqrt(bigfac(n)) * sqrt(bigfac(m))))*
    sum([bigbinomial(n,k)*bigbinomial(m,k) * (1/(bigpw2(k))) * bigfac(k) * (im*big(ξ0)/2)^(n + m - 2*k) for k = 0:minimum([n,m])])

Ty(n::Int64, m::Int64, ξ0::Float64) = exp(- ξ0^2 / 4) * sqrt(bigfac(n))/sqrt(bigfac(m)) * big(-ξ0/sqrt(2))^(m-n) * mylaguerre2big(m-n, n, big(ξ0^2 /2))

# work in basis (n,m,ky0) that goes as (n=0,ky=ky0), (0,ky0+2π/a), ..., (0,ky0+(p-1)2π/a), (1,ky0),.....
# hamiltonian is composed of A(on diagional) and B matrices that are nonzero on the 3 diagonals

# block matrices on the diagonal -- correspond to (n,k)*(n,k') elements
function matA(n::Int64, ξ0::Float64, ky0::Float64, Y::Float64, U0::Float64, a::Float64, p::Int64)
    # get elements from above for the n-th LL
    Tnn = Ty(n, n, ξ0)
    En = E_LL(n, ξ0, a)

    # diagonal elements
    d(j) = En + U0 * Tnn * cos(ξ0^2 * a /(2π) * (ky0 + j * 2π/a))

    # elements on the two second diagonals
    f = exp(-im*Y/p) * U0 * Tnn / 2

    #construct matrix
    mat::Matrix{ComplexF64} = diagm(0 => [d(j) for j = 0:(p-1)], 1 => f * ones(p-1), -1 => conj(f) * ones(p-1))
    mat[1, p] += conj(f)
    mat[p, 1] += f

    return mat
end

# matrix elements away from diagonal; LL mixing
function matB(n::Int64, m::Int64, ξ0::Float64, ky0::Float64, Y::Float64, U0::Float64, a::Float64, p::Int64)
    # get elements from above for LL n and m
    Txnm = Tx(n, m, ξ0)
    Txmn_conj = conj(Txnm)
    Tynm = Ty(n, m, ξ0)
    Tymn = Ty(m, n, ξ0)

    # diagonal elements
    d(j) = U0/2 * (Txnm * exp(im * ξ0^2 * a /(2π) * (ky0 + j * 2π/a)) + Txmn_conj * exp(-im * ξ0^2 * a /(2π) * (ky0 + j * 2π/a)) )

    # elements on the two second diagonals
    f1 = exp(-im*Y/p) * U0 * Tynm / 2
    f2 = exp(im*Y/p) * U0 * Tymn / 2

    #construct matrix
    mat::Matrix{ComplexF64} = diagm(0 => [d(j) for j = 0:(p-1)], 1 => f1 * ones(p-1), -1 => f2 * ones(p-1))
    mat[1, p] += f2
    mat[p, 1] += f1

    return mat
end

# row n of the hamiltonian matrix
function get_ham_row(n::Int64, ξ0::Float64, ky0::Float64, Y::Float64, U0::Float64, a::Float64, p::Int64, NLL::Int64)
    if n == 0
        arrayB2 = reduce(hcat, [matB(n, m, ξ0, ky0, Y, U0, a, p) for m = n+1:NLL])
        return [matA(n, ξ0, ky0, Y, U0, a, p) arrayB2]
    elseif n == NLL
        arrayB1 = reduce(hcat, [matB(n, m, ξ0, ky0, Y, U0, a, p) for m = 0:n-1])
        return [arrayB1 matA(n, ξ0, ky0, Y, U0, a, p)]
    else
        arrayB1 = reduce(hcat, [matB(n, m, ξ0, ky0, Y, U0, a, p) for m = 0:n-1])
        arrayB2 = reduce(hcat, [matB(n, m, ξ0, ky0, Y, U0, a, p) for m = n+1:NLL])
        return [arrayB1 matA(n, ξ0, ky0, Y, U0, a, p) arrayB2]
    end    
end

# full hamiltonian matrix
function get_full_ham(ξ0::Float64, ky0::Float64, Y::Float64, U0::Float64, a::Float64, p::Int64, NLL::Int64)
    rows_vector = [get_ham_row(n, ξ0, ky0, Y, U0, a, p, NLL) for n = 0:NLL]
    ham = reduce(vcat, rows_vector)
    return Hermitian(ham)
end



end