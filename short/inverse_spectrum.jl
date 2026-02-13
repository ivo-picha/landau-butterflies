using LinearAlgebra
using Plots
using Primes
using ProgressMeter
using Measures

# parameters
a = 50e-10; # lattice constant [m]
a_in_angstr = a*1e10
U0 = 0.01; # strength of periodic potential [eV]

# physical constants
ħ = 6.62607015e-34/(2π); # Planck constant [J s]
e = 1.602176634e-19; # elementary charge [C]
m_e = 9.1093837139e-31; # electron mass [kg];

#number of bands to consider 
p_set = 301
if !isprime(p_set)
    println("would be nice if you chose a prime p")
end

# list of q values to iterate over
Nq = 400; #number of desired steps of q 
q_list = unique(Int.(round.(range(p_set, 3*p_set, Nq))))

#set up list of ky* values to iterate over
Nky = 65;
ky_list = range(0, 2π/a, Nky)

# empty lists to store coordinates for plot
energies = [];
phis = [];

@showprogress for qn in q_list
    p = p_set
    q = qn
    # gcdpq = gcd(p,q)
    # if gcdpq != 1 && p == q # make sure only coprime fractions enter
    #     p = Integer(p / gcdpq)
    #     q = Integer(q / gcdpq)
    # end

    #magnetic field dependent parameters
    phi = p/q # unit flux per unit cell
    phi_inv = q/p
    xi0 = sqrt(2π / phi)
    lB = xi0 * a / (2π) # magnetic length
    E0 = ħ^2 / (2 * e * m_e * lB^2) #lowest landau level energy
    T0 = exp(-xi0^2 / 4)

    B = U0 * T0 / 2 # non-diagonal matrix element

    # println(E0)
    # println(B)

    xi0 = sqrt(2π / phi)
    lB = xi0 * a / (2π) # magnetic length
    E0 = ħ^2 / (2 * e * m_e * lB^2)

    for ky in ky_list

        A(m) = E0 + U0 * T0 * cos(xi0 * lB * (ky + m * 2π / a)) # diagonal matrix element

        #construct Hamiltonian
        H = diagm(0 => [A(m) for m = 1:p], 1 => B * ones(p-1), -1 => B * ones(p-1))
        H[1,p] = B
        H[p,1] = B

        H = Hermitian(H)

        evalsH = eigvals(H) .- E0 # remove cyclotron frequency offset


        energies = [energies; evalsH] #add eigenvalues to list of energies

    end

    phis = [phis; phi_inv * ones(p * length(ky_list))] #add phi values to list of phis

end


Npoints = length(energies)
println("plotting $Npoints number of points")

plot1 = scatter(phis, energies,title = "U₀ = $U0 eV, a = $a_in_angstr Å", markersize = 0.4, color = :black, label = "", 
    xlabel = "1/ϕ = q/p", ylabel = "E - ε(ϕ) [eV]", framestyle = :box, size = (1000,500), tickfontsize = 16, guidefontsize = 18, margin=9mm,
    yticks = false);


# plot enveloping function the lazy way
# function func1(pi1)
#     phi = 1/pi1
#     xi0 = sqrt(2π / phi)
#     lB = xi0 * a / (2π) # magnetic length
#     E0 = ħ^2 / (2 * e * m_e * lB^2)
#     T0 = exp(-xi0^2 / 4)
#     return (E0 + 2 *U0 * T0) 
# end

# pi1_grid = range(1.,2.,200)
# envss = func1.(pi1_grid)

# plot!(plot1, pi1_grid, envss, color = :red, lw = 1.5, label = "ϵ₀ + 2 U₀ Θ", legendfontsize = 16);

savefig(plot1, "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/inverse_spectrum.png")