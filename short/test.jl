using Plots

q = [1,2,3,4,5,6].+0.5
Uc = [0.0151,0.0083,0.0058,0.0047,0.0040,0.0035]


scatter(log.(q),log.(Uc); 
    xlabel="q", ylabel="Uc*q", title="Critical potential", markershape=:diamond, label="vLL", color=:red)

x = log.(q)
y = log.(Uc)
X = hcat(ones(length(x)), x)
coef = X \ y
intercept, slope = coef[1], coef[2]
plot!(x, X * coef, label="slope= $slope", color=:red, lw=2)



