using Plots

# q = [1,2,3,4,5]
# Uc = [0.0151,0.0083,0.00585,0.0047,0.0040]

# scatter([1/qi for qi in q], Uc; 
#     xlabel="1/q", ylabel="Uc", title="Critical potential", markershape=:diamond, label="vLL", color=:red, xlims=(0,1), ylims=(0,0.02))


# scatter(log.(q),log.(Uc); 
#     xlabel="ln q", ylabel="ln Uc", title="Critical potential", markershape=:diamond, label="vLL", color=:red)

# x = log.(q)
# y = log.(Uc)
# X = hcat(ones(length(x)), x)
# coef = X \ y
# intercept, slope = coef[1], coef[2]
# plot!(x, X * coef, label="slope= $(round(slope;digits=4))", color=:red, lw=2)

# plot!(x,X*[coef[1],-1.0])


# include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
# using .Hamil 


[Hamil.T(q-1,q-1,Float32.(1/q))-Hamil.T(q,q,Float32.(1/q)) for q in 1:10]

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

# matrix elements (Θ in overleaf)
#T(n::Int64, m::Int64, phi::Float32) = Float32(exp(-Float32(π)/(Float32(2.0)*phi)) * sqrt(bigfac(min(n,m))/bigfac(max(n,m))) * (Float32(π)/phi)^(abs(n-m)/2) * mylaguerre2big(abs(n-m),min(n,m),Float32(π)/phi))
function T(n::Int64, m::Int64, phi::Float32)
    tt = exp(-Float32(π)/(Float32(2.0)*phi))
    tt *= sqrt_factorial_ratio(n,m)
    tt *= laguerre_scaled(n, m, Float32(π)/phi)
    return Float32(tt)
end

T(20,300,1f0)