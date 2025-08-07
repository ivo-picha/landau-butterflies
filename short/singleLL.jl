using LinearAlgebra
using Plots
using Primes
using ProgressMeter
using Measures
using Polynomials, SpecialPolynomials

# parameters
a = 5e-9; # lattice constant [m]
U0 = 0.05; # strength of periodic potential [eV]

# physical constants
ħ = 6.62607015e-34/(2π); # Planck constant [J s]
e = 1.602176634e-19; # elementary charge [C]
m_e = 9.1093837139e-31; # electron mass [kg];

# n-th Landau level (0 lowest)
nLL = 0;

#number of bands to consider 
q = 120

# list of q values to iterate over
startphi = 1.0;
endphi = 4.0;
p_list = Int.(collect(range(round(q*startphi),round(q*endphi))))

#set up list of ky* values to iterate over
Nky = 10;
ky_list = range(0, 2π/a, Nky)
ky_list = ky_list[1:end-1] # 0 = 2π ; remove repetition

NY = Nky;
Y_list = range(0, 2π, NY)
Y_list = Y_list[1:end-1] # 0 = 2π ; remove repetition

# define an easier to use laguerre polynomial 
function mylaguerre(α::Number, n::Integer, x::Number)
    list1 = push!(zeros(Integer,n),1)
    lag = Laguerre{α}(list1)
    return lag(x)    
end

# LL energy in eV
E_LL(ξ0::Float64) = (nLL + 0.5) * ħ^2 / (e * m_e * (ξ0 * a / (2π))^2)

#matrix elements
T(ξ0::Float64) = exp(- ξ0^2 / 4) * mylaguerre(0, nLL, ξ0^2 /2)

# empty lists to store coordinates for plot
energies = [];
phis = [];

@showprogress for p in p_list

    pj = Int(p/gcd(p,q))
    qj = Int(q/gcd(p,q))
    #magnetic field dependent parameters
    phi = pj/qj # unit flux per unit cell
    xi0 = sqrt(2π / phi)
    lB = xi0 * a / (2π) # magnetic length
    En = E_LL(xi0) #landau level energy
    Tn = T(xi0)

    B = U0 * Tn / 2 # non-diagonal matrix element

    for ky in ky_list
        for Y in Y_list

        A(m) = En + U0 * Tn * cos(xi0 * lB * (ky + m * 2π / a)) # diagonal matrix element

        #construct Hamiltonian
        H = diagm(0 => [A(m) for m = 1:pj], 1 => B* exp(-im*Y/pj) * ones(pj-1), -1 => B * exp(im*Y/pj) * ones(pj-1))
        H[1,pj] += B * exp(im*Y/pj)
        H[pj,1] += B * exp(-im*Y/pj)

        H = Hermitian(H)

        evalsH = eigvals(H)

        energies = [energies; evalsH] #add eigenvalues to list of energies
        phis = [phis; phi * ones(length(evalsH))] #add phi values to list of phis

        end
    end
end

Npoints = length(energies)
println("plotting $Npoints number of points")

plot1 = scatter(phis, energies,title = "U₀ = $U0 eV, a = 50 Å, n = $nLL", markersize = 0.4, color = :black, label = "", 
    xlabel = "ϕ = p/q", ylabel = "E [eV]", framestyle = :box, size = (1200,800), tickfontsize = 16, guidefontsize = 18, margin=9mm)


# plot enveloping functions
# phi_list_unique = unique(phis)
# xi0list = sqrt(2π)./ sqrt.(phi_list_unique)
# y_LLenergies = E_LL.(xi0_list_denser)

# plot!(plot1, phi_inv_list_denser, y_env_upper, color = :red, lw = 1.5, linestyle = :dash, label = "", legendfontsize = 16);
# plot!(plot1, phi_inv_list_denser, y_env_lower, color = :red, lw = 1.5, linestyle = :dash, label = "ϵₙ ± 2 U₀ Θₙₙ", legendfontsize = 16);
# plot!(plot1, phi_inv_list_denser, y_LLenergies, color = :red, lw = 1.5, label = "ϵₙ", legendfontsize = 16);

savefig(plot1, "plots/local/singleLLs/n$nLL-pi$startphi-pf$endphi-U$U0.png")