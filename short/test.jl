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


include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil 
using LinearAlgebra
eigvals(Hamil.get_full_ham(0.1f0, 0f0, 0f0, 0.01f0, 5f-9, 1, 1500))



n(idx,p) = fld(idx, p)
m(idx,p) = idx%p

pp = 5
for idx in 1:(pp*3)
    println("idx=$idx -> n=$(n(idx,pp)), m=$(m(idx,pp))")
end