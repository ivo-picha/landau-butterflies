using Plots

#list of nonrepeating p and q that were calculated
# list_p = [1,2,3,4,5,6,7,8]
# list_q = copy(list_p)
# ip = reshape(collect(Base.product(list_p, list_q)),:)
# past_phi = []
# pqs = Tuple{Int64,Int64}[]
# for (pj,qj) in ip
#     phi = Float64(pj/qj)

#     if in(phi,past_phi)
#         continue
#     elseif phi > 3
#         continue
#     end

#     push!(past_phi,phi)
#     push!(pqs,(pj,qj))

# end
# mask = sortperm(past_phi)
# pqs_ord = pqs[mask]
# println(pqs_ord)
(0.000273+0.000241+0.000262+0.000341+0.000301)/5

pqlist = [(0,1),
          (1, 8), (1, 7), (1, 6), (1, 5), (1, 4), (2, 7), (1, 3), (3, 8), (2, 5), (3, 7), 
          (1, 2), (4, 7), (3, 5), (5, 8), (2, 3), (5, 7), (3, 4), (4, 5), (5, 6), (6, 7),
          (7, 8), (1, 1), (8, 7), (7, 6), (6, 5), (5, 4), (4, 3), (7, 5), (3, 2), (8, 5),
          (5, 3), (7, 4), (2, 1), (7, 3), (5, 2), (8, 3), (3, 1)]

# all values in lists in eV; t's with minus signs
#--------------- U0 = 20 meV
# chemical potential
mu2 = [0.00026486,
       missing, 0.0002836, 0.000284, 0.000304, 0.000313, missing, 0.000359, missing, 0.0004156, missing,
       0.000515, missing, 0.000645, missing, 0.000745, missing, 0.000878, 0.000961, 0.00101, missing,
       missing, 0.001319, missing, 0.00176, missing, missing, 0.002339, missing, 0.003001, 0.003416,
       0.003702, 0.00407, 0.005245, 0.00699, 0.007916, 0.008877, 0.010872]
# NN hopping
tnn2 = [0.00183264,
        missing, 0.001787, 0.00178, 0.001763, 0.001744, missing, 0.001709, missing, 0.001652, missing,
        0.001632, missing, 0.001592, missing, 0.001557, missing, 0.001497, 0.001442, 0.001402, missing,
        missing, 0.001178, missing, 0.000974, missing, missing, 0.000866, missing, 0.000785, 0.000761,
        0.000723, 0.000683, 0.00055, 0.000407, 0.000354, 0.000308, 0.000225]
# NNN hopping
tnnn2 = [0.00022359,
         missing, missing, 0.000205, 0.000203, 0.000203, missing, 0.000201, missing, 0.0002, missing,
         0.00017, missing, 0.00014, missing, 0.000145, missing, 0.000172, 0.000143, 0.000175, missing,
         missing, 0.000123, missing, 6.1e-5, missing, missing, 6.7e-5, missing, 5.8e-5, 5.2e-5,
         4.9e-5, 4.6e-5, 3.3e-5, 1.9e-5, 1.6e-5, 1.2e-5, 7.0e-6] 
# d-NNN hopping
tdnnn2 = [1.4310611595233178e-5,
          missing, missing, 6.0e-5, 9.2e-5, 0.000128, missing, 0.000194, missing, 0.000531, missing,
          0.000323, missing, 0.000454, missing, 0.000979, missing, 0.000382, 0.000318, 0.000241, missing,
          missing, 0.000211, missing, 6.6e-5, missing, missing, 9.1e-5, missing, 0.000167, 0.000106,
          0.000124, 0.0001, 6.5e-5, 1.2e-5, 1.75e-5, 2.2e-5, 1.7e-5]

#--------------- U0 = 50 meV
# chemical potential
mu5 = [-0.031030036,
       missing, missing, -0.030994, -0.031, -0.030989, missing, -0.030941, missing, -0.030909, missing,
       -0.030808, missing, -0.030757, missing, -0.03069, missing, -0.030598, missing, -0.030588, missing,
       missing, -0.030194, missing, -0.029987, missing, missing, -0.029655, missing, -0.029295, -0.029079,
       -0.02891, -0.0287, -0.027963, -0.02695, -0.026367 , -0.025772, -0.024429]
# NN hopping
tnn5 = [0.00061797,
        missing, missing, 0.000645, 0.00064, 0.000638, missing, 0.000632, missing, 0.000625, missing,
        0.000615, missing, 0.000602, missing, 0.000592, missing, 0.000582, missing, 0.000567, missing,
        missing, 0.000537, missing, 0.000497, missing, missing, 0.00046, missing, 0.000424, 0.000405,
        0.00039, 0.000373, 0.000315, 0.00026, 0.000233, 0.000206, 0.000157]
# NNN hopping
tnnn5 = [1.5045e-5,
         missing, missing, 2.0e-5, 1.4e-5, 2.0e-5, missing, 2.1e-5, missing, 2.0e-5, missing,
         1.7e-5, missing, 1.7e-5, missing, 1.4e-5, missing, 1.65e-5, missing, 1.7e-5, missing,
         missing, 2.2e-5, missing, 1.3e-5, missing, missing, 1.1e-5, missing, 3.0e-6, 9.0e-6,
         9.0e-6, 9.0e-6, 6.0e-6, 2.0e-6, 5.0e-6, 3.0e-6, 8.0e-6]
# d-NNN hopping
tdnnn5 = [1.473793319102998e-6,
          missing, missing, 9.0e-6, 1.5e-5, 9.0e-6, missing, 1.4e-5, missing, 0.000152, missing,
          2.6e-5, missing, 5.8e-5, missing, 0.00019, missing, 4.1e-5, missing, 7.0e-5, missing,
          missing, 2.5e-5, missing, 1.2e-5, missing, missing, 3.0e-5, missing, 0.000176, 1.7e-5,
          2.0e-5, 2.0e-5, 1.3e-5, 2.0e-6, 2.0e-6, 6.0e-6, 2.0e-6] 

#--------------- U0 = 80 meV
# chemical potential
mu8 = [-0.0701928,
       missing, missing, -0.07004, -0.0703, -0.07014, missing, -0.07011, missing, -0.070094, missing,
       -0.070028, missing, -0.069975, missing, -0.069925, missing, -0.06986, missing, missing, missing,
       missing, -0.069422, missing, missing, -0.06936, missing, -0.069156, missing, -0.068867, missing,
       missing, missing, -0.06779, missing, -0.066672, missing, -0.065098]
# NN hopping
tnn8 = [0.00016716,
        missing, missing, missing, 0.000241, 0.000248, missing, 0.000246, missing, 0.000244, missing,
        0.000241, missing, 0.000236, missing, 0.000235, missing, 0.000234, missing, missing, missing,
        missing, 0.000201, missing, missing, missing, missing, 0.000194, missing, 0.000185, missing,
        missing, missing, 0.000164, missing, 0.000113, missing, 7.9e-5]
# NNN hopping
tnnn8 = [7.949579256338586e-5,
         missing, missing, missing, 4.0e-6, 5.0e-6, missing, 2.0e-6, missing, 2.0e-6, missing,
         2.0e-6, missing, 8.0e-6, missing, 2.0e-6, missing, 5.0e-6, missing, missing, missing,
         missing, 5.2e-5, missing, missing, missing, missing, 4.0e-6, missing, 7.0e-6, missing,
         missing, missing, 1.7e-5, missing, 1.6e-5, missing, 2.3e-5] 
# d-NNN hopping
tdnnn8 = [4.814433585517812e-7,
          missing, missing, missing, 1.7e-5, 1.0e-6, missing, 2.0e-6, missing, 2.8e-5, missing,
          0.0, missing, 1.2e-5, missing, 9.0e-5, missing, 6.0e-6, missing, missing, missing,
          missing, 1.4e-5, missing, missing, missing, missing, 1.3e-5, missing, 2.0e-5, missing,
          missing, missing, 6.0e-6, missing, 1.0e-6, missing, 2.7e-5] 


## y NN hopping phases, multiples of pi; for U = 50 meV only
ph_list = [
    # q = 1
    (1, 1, [1-0.012694]),
    (2, 1, [-0.990117]),
    (3, 1, [1-0.002523]),

    # q = 2
    (1, 2, [0.496945, -0.496481]),
    (3, 2, [-0.499666, 0.499569]),
    (5, 2, [-0.497311, 0.497408]),

    # q = 3
    (1, 3, [-0.667798, 0.667652, -3.7e-5]),
    (4, 3, [0.999982, -0.331544, 0.331336]),
    (5, 3, [0.667063, -0.66709, -5.6e-5]),
    (7, 3, [-0.665525, 0.668157, -0.002667]),
    (8, 3, [0.999956, -0.333524, 0.333499]),

    # q = 4
    (1, 4, [-0.748665, 0.748629, -0.249971, 0.250004]),
    (3, 4, [-0.75159, 0.751568, -0.251694, 0.251126]),
    (5, 4, [missing, missing, missing, missing]),
    (7, 4, [-0.750491, 0.750438, -0.250012, 0.249871]),

    # q = 5
    (1, 5, [-0.804509, 0.797776, 0.204444, -0.401694, 0.203342]),
    (3, 5, [0.800265, -0.800088, 0.400162, -0.398982, 4.9e-5]),
    (4, 5, [missing, missing, missing, missing, missing]),
    (6, 5, [missing, missing, missing, missing, missing]),
    (7, 5, [-0.819819, 0.783955, 0.387833, 0.338673, -0.001025]),
    (8, 5, [-0.999983, -0.600326, 0.600311, -0.200033, 0.200124]),

    # q = 6; all sus
#     (1, 6, [-0.833335, 0.666772, -0.499989, 0.665879, -0.168557, 0.166653]),
#     (5, 6, [-0.765357, 0.825509, -0.299166, 0.486481, 0.168008, -0.324662]),
#     (7, 6, [-0.833287, 0.833438, -0.499896, 0.331251, 0.333338, -0.165692]),

    # q = 7 ; all sus
#     (1, 7, fill(missing, 7)),
#     (2, 7, fill(missing, 7)),
#     (3, 7, fill(missing, 7)),
#     (4, 7, fill(missing, 7)),
#     (5, 7, [0.855682, -0.856724, -0.557214, 0.560345, -0.201497, -0.195273, -6.4e-5]),
#     (6, 7, [0.85814, -0.717444, 0.855484, -0.130726, -0.426493, 0.430529, 0.142949]),
#     (8, 7, fill(missing, 7)),
]

θ = range(0, 2π, length=500)
circles = plot(
    size = (350,300),
    framestyle = :box,
    xticks=false,
    yticks=false,
    legend=:outerright,
)
for n in 1:5
    r = 0.01 + n / 4
    plot!(
        circles,
        r .* cos.(θ),
        r .* sin.(θ),
        label="",
        color = :black
    )
end

function plot_phases!(plot::Plots.Plot, set::Tuple{Int64,Int64,Vector})
        p, q, pts = set
        r = 0.01 + q / 4
        ptsxy = [(r*cos(pt*π), r*sin(pt*π)) for pt in pts]
        scatter!(plot, ptsxy, color = (iseven(p) ? :red : :blue), label = "")
end

for set in ph_list
        plot_phases!(circles, set)
end
savefig




## ---------------------- plotting 

pqplotlist = [(1, 8), (1, 3), (1, 2), (2, 3),(5, 6), (1, 1), (7,6), (4, 3),  (3, 2), (7, 4), (2, 1), (7, 3), (5, 2), (8, 3), (3, 1)]
phiplotlist = [pj/qj for (pj,qj) in pqplotlist]
labellist = ["$pj/$qj" for (pj,qj) in pqplotlist]
pushfirst!(phiplotlist, 0)
pushfirst!(labellist, "0")

philist = [pj/qj for (pj,qj) in pqlist]

plt_opt = (
    framestyle = :box,
    xaxis = "ϕ = p/q",
    xlims = (0,3),
    size = (480,320),
    xticks = (phiplotlist,labellist)
)

plot_mu = scatter(philist,mu2.*1000, label= "U₀ = 20 meV";plt_opt...)
scatter!(plot_mu,philist,mu5.*1000, label= "U₀ = 50 meV";plt_opt...)
scatter!(plot_mu,philist,mu8.*1000, label= "U₀ = 80 meV";plt_opt...)
yaxis!(plot_mu, "μ [meV]")

plot_tnn = scatter(philist,tnn2.*1000, label= "";plt_opt...)
scatter!(plot_tnn,philist,tnn5.*1000, label= "";plt_opt...)
scatter!(plot_tnn,philist,tnn8.*1000, label= "";plt_opt...)
yaxis!(plot_tnn, "|t| [meV]")

plot_tnnn = scatter(philist,tnnn2.*1000, label= "";plt_opt...)
scatter!(plot_tnnn,philist,tnnn5.*1000, label= "";plt_opt...)
scatter!(plot_tnnn,philist,tnnn8.*1000, label= "";plt_opt...)
yaxis!(plot_tnnn, "|t'| [meV]")

plot_tdnnn = scatter(philist,tdnnn2.*1000, label= "";plt_opt...)
scatter!(plot_tdnnn,philist,tdnnn5.*1000, label= "";plt_opt...)
scatter!(plot_tdnnn,philist,tdnnn8.*1000, label= "";plt_opt...)
yaxis!(plot_tdnnn, "|t''| [meV]")

plt_combi = plot(plot_mu, plot_tnn, plot_tnnn, plot_tdnnn, layout = (2,2), size = (1.6*600,500))