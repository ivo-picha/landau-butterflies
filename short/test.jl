# define a function which can handle big factorials beyond Int128
function bigfac(n::Int)
    if n > 20
        return factorial(big(n))
    else
        return factorial(n)
    end    
end

function bigbinomial(n::Int, m::Int)
        if n > 60
        return binomial(big(n),big(m))
    else
        return binomial(n,m)
    end    
end

function bigpw2(n::Int)
        if n > 60
        return big(2)^n
    else
        return 2^n
    end    
end

Tx(n::Int64, m::Int64, ξ0::Float64) = exp(- ξ0^2 / 4) * (sqrt((bigpw2(n+m)))/(sqrt(bigfac(n)) * sqrt(bigfac(m))))*
    sum([bigbinomial(n,k)*bigbinomial(m,k) * (1/(bigpw2(k))) * bigfac(k) * (im*big(ξ0)/2)^(n + m - 2*k) for k = 0:minimum([n,m])])

for n=0:100
    for m=0:100
        # if n<m
        #     continue
        # end
        outp = Tx(n,m,0.01)
        println("working on n=$n, m=$m; Tx outputs $outp")
    end
end