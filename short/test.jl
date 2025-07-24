using Polynomials, SpecialPolynomials

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


# matrix elements (Θ in overleaf)
Tx(n::Int64, m::Int64, ξ0::Float64) = exp(- ξ0^2 / 4) * (sqrt((bigpw2(n+m)))/(sqrt(bigfac(n)) * sqrt(bigfac(m))))*
    sum([bigbinomial(n,k)*bigbinomial(m,k) * (1/(bigpw2(k))) * bigfac(k) * (im*big(ξ0)/2)^(n + m - 2*k) for k = 0:minimum([n,m])])

Ty(n::Int64, m::Int64, ξ0::Float64) = exp(- ξ0^2 / 4) * sqrt(bigfac(n))/sqrt(bigfac(m)) * big(-ξ0/sqrt(2))^(m-n) * mylaguerre2big(m-n, n, big(ξ0^2 /2))


x0=sqrt(2π/0.3)
limn = 1e5
txn = NaN
tyn = NaN
g=100
for n=1:g
    for m = 1:g
        txn = Tx(n,m,x0)
        tyn = Ty(n,m,x0) # starts overflowing later
        if abs(txn) > limn || abs(tyn) > limn
            println("overflow at n=$n, m=$m, giving Tx = $txn, Ty = $tyn")
        end
    end
end



# nn = 30
# mm = 2


# #exp(- x0^2 / 4) * sqrt(bigfac(nn))/sqrt(bigfac(mm)) * big(x0/sqrt(2))^(mm-nn) * (-1.0)^(nn-mm) * mylaguerre2big(mm-nn, nn, x0^2 /2)

# (mylaguerre2big(mm-nn, nn, x0^2 /2), mylaguerre(mm-nn, nn, x0^2 /2))







# (im*big(x0)/2)^(nn + mm)
# s = sum([bigbinomial(nn,k)*bigbinomial(mm,k) * (1/(bigpw2(k))) * bigfac(k) * (im*big(x0)/2)^(nn + mm - 2*k) for k = 0:minimum([nn,mm])])

# bs = exp(- x0^2 / 4) * (sqrt((bigpw2(nn+mm)))/(sqrt(bigfac(nn)) * sqrt(bigfac(mm))))

# # exp(- x0^2 / 4) * (sqrt((bigpw2(nn+mm)))/(sqrt(bigfac(nn)) * sqrt(bigfac(mm))))*
# #     sum([bigbinomial(nn,k)*bigbinomial(mm,k) * (1/(bigpw2(k))) * bigfac(k) * (im*big(x0)/2)^(nn + mm - 2*k) for k = 0:minimum([nn,mm])])

# println("sum is $s")
# println("prefactor is $bs")
# println("Tx is $(s*bs)")