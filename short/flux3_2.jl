using LinearAlgebra
using Plots


function hamil32(X,Y)
    h = diagm([2*cos(X/3 + 4π/3*m) for m=0:2])
    h += diagm(1=>[exp(-im*Y/3) for m=0:1])
    h[1,3] += exp(im*Y/3)
    return Hermitian(h)
end

XYrange = [(t*2π,0) for t=0:0.01:1]
XYrange = [XYrange; [(2π,t*π) for t=0:0.01:0.99]]
append!(XYrange, [(2π*(1-t),π*(1-t)) for t=0:0.01:0.99])

energies = Float64[];
js = Float64[];
for (j,(X,Y)) in enumerate(XYrange)
    append!(energies,eigvals(hamil32(X,Y)))
    append!(js,j)
end

scatter(js,energies)