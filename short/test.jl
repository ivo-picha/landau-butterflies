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

Float32(π)^(70)