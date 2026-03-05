using Plots
using LinearAlgebra
using Measures

include(joinpath(dirname(@__DIR__),"funcs/hamiltonian.jl"))
using .Hamil                        # build a Hamiltonian matrix in a Landau level basis

q = collect(1:50)
Uc = [0.0151,0.0083,0.0059,0.0047,0.0040,0.0035,0.0033,0.0032,0.0031,0.0033,0.0035,0.0036,0.0038,0.0042,0.0048,0.0052,0.0059,0.0064,0.0066,0.0069,0.0072,0.0075,0.0076,0.0078,0.008,0.0082,0.0083,0.0084,0.0085,0.0085,0.0086,0.0087,0.0088,0.0089,0.0089,0.0089,0.009,0.0091,0.0091,0.0091,0.0092,0.0092,0.0093,0.0093,0.0093,0.0093,0.0094,0.0094,0.0095,0.0095]

function get_Uc_2b(q)
    phi = Float32(1/q)
    a = 5f-9
    Uc = Hamil.E_LL(q,phi,a) - Hamil.E_LL(q-1,phi,a)
    Uc /= 2f0
    Uc /= Hamil.T(q,q,phi) - Hamil.T(q-1,q-1,phi)
    return abs(Uc)
end

Uc_2b = get_Uc_2b.(q)

# scatter([1/qi for qi in q], Uc; 
#     xlabel="1/q", ylabel="Uc", title="Critical potential", markershape=:diamond, label="vLL", color=:red, xlims=(0,1), ylims=(0,0.02))
tf = 18
lf = 20
plt = plot(xlabel="q", ylabel="log(Critical U₀)", framestyle=:box,size = (950,500), legend=:topleft, xtickfont = tf, ytickfont = tf, guidefont = lf, legendfont = tf-5, margin = 6mm)

plot!(plt,q,log.(Uc_2b);
    color=1, label="2-LL; analytical", markershape=:diamond, linestyle = :dash, lw=2, ms = 8)

plot!(plt,q,log.(Uc); 
    markershape=:diamond, label="all-LL; numerical", color=2, lw=2, ms = 8)

# scatter(log.(q),log.(Uc); 
#     xlabel="ln q", ylabel="ln Uc", title="Critical potential", markershape=:diamond, label="vLL", color=:red)

# zero field critical potential
Uczf = 0.0106
hline!(plt, [log(Uczf)], label="ϕ=0", color=:black, linestyle=:dot, lw=2)


## 
# function f4(q::Integer,U0::Real)
#     if q<2
#         error("q must be at least 2 for 4 band calculation")
#     end
#     phi = Float32(1/q)
#     a = 5f-9
#     pm = (-1)^(q)
#     e1 = Hamil.E_LL(q-2,phi,a) 
#     e2 = Hamil.E_LL(q-1,phi,a)
#     e3 = Hamil.E_LL(q,phi,a) 
#     e4 = Hamil.E_LL(q+1,phi,a)
#     #
#     t11 = Hamil.T(q-2,q-2,phi)
#     t22 = Hamil.T(q-1,q-1,phi)
#     t33 = Hamil.T(q,q,phi)
#     t44 = Hamil.T(q+1,q+1,phi)
#     #
#     t13 = Hamil.T(q,q-2,phi)
#     t24 = Hamil.T(q-1,q+1,phi)
#     f = -e1+e2-e3+e4
#     f += -sqrt(e1^2 + e3^2 - 2*e1*(e3 + pm*2*U0*(t11-t33)) + 4*U0^2 *(4*t13^2 + (t11-t33)^2) + pm*4*U0*e3*(t11-t33))
#     f += -sqrt(e2^2 + e4^2 - 2*e2*(e4 + pm*2*U0*(t22-t44)) + 4*U0^2 *(4*t24^2 + (t22-t44)^2) + pm*4*U0*e4*(t22-t44))
#     f += pm*2*U0*(t11-t22+t33-t44)
#     return f/2.0
# end

# tr = range(0.0,0.016,1000)
# fr = f4.(3,tr)
# plot(tr,fr)



# function get_ev(q::Integer,V::Real)
#     phi = Float32(1/q)
#     a = 5f-9

#     qrange = collect(q-2:q+1)

#     M = reshape(round.([cos(π*abs(n-m)/2) for n in qrange for m in qrange]; digits =3),4,4)
#     T = Array{Float32}(undef,4,4)
#     for (n,qn) in enumerate(qrange)
#         for (m,qm) in enumerate(qrange)
#             T[n,m] = Hamil.T(qn,qm,phi)
#         end
#     end
#     E = diagm([Hamil.E_LL(qn,phi,a) for qn in qrange])
#     H = E .+ (2.0*V*(-1)^(q+1)).* (M .* T)
#     H = Hermitian(H)
#     return eigvals(H)
# end

# Vlist = range(0f0,0.016f0,500)
# elist = get_ev.(3,Vlist)
# ehlist = hcat(elist...)'
# plot(Vlist,ehlist)



