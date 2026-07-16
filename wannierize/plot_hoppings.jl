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
(0.000943+0.000979+0.000127+0.000961+0.001058+0.001983)/6

pqlist = [(1, 8), (1, 7), (1, 6), (1, 5), (1, 4), (2, 7), (1, 3), (3, 8), (2, 5), (3, 7), 
          (1, 2), (4, 7), (3, 5), (5, 8), (2, 3), (5, 7), (3, 4), (4, 5), (5, 6), (6, 7),
          (7, 8), (1, 1), (8, 7), (7, 6), (6, 5), (5, 4), (4, 3), (7, 5), (3, 2), (8, 5),
          (5, 3), (7, 4), (2, 1), (7, 3), (5, 2), (8, 3), (3, 1)]

# all values in lists in eV; t's with minus signs
#--------------- U0 = 20 meV
# chemical potential
mu2 = [missing, missing, 0.000284, 0.000304, 0.000313, missing, 0.000359, missing, 0.0004156, missing,
       0.000515, missing, 0.000645, missing, 0.000745, missing, 0.000878, 0.000961, 0.00101, missing,
       missing, 0.001319, missing, 0.00176, missing, missing, 0.002339, missing, 0.003001, 0.003416,
       0.003702, 0.00407, 0.005245, 0.00699, 0.007916, 0.008877, 0.010872]
# NN hopping
tnn2 = [missing, missing, 0.00178, 0.001763, 0.001744, missing, 0.001709, missing, 0.001652, missing,
        0.001632, missing, 0.001592, missing, 0.001557, missing, 0.001497, 0.001442, 0.001402, missing,
        missing, 0.001178, missing, 0.000974, missing, missing, 0.000866, missing, 0.000785, 0.000761,
        0.000723, 0.000683, 0.00055, 0.000407, 0.000354, 0.000308, 0.000225]
# NNN hopping
tnnn2 = [missing, missing, 0.000205, 0.000203, 0.000203, missing, 0.000201, missing, 0.0002, missing,
         0.00017, missing, 0.00014, missing, 0.000145, missing, 0.000172, 0.000143, 0.000175, missing,
         missing, 0.000123, missing, 6.1e-5, missing, missing, 6.7e-5, missing, 5.8e-5, 5.2e-5,
         4.9e-5, 4.6e-5, 3.3e-5, 1.9e-5, 1.6e-5, 1.2e-5, 7.0e-6] 
# d-NNN hopping
tdnnn2 = [missing, missing, 6.0e-5, 9.2e-5, 0.000128, missing, 0.000194, missing, 0.000531, missing,
          0.000323, missing, 0.000454, missing, 0.000979, missing, 0.000382, 0.000318, 0.000241, missing,
          missing, 0.000211, missing, 6.6e-5, missing, missing, 9.1e-5, missing, 0.000167, 0.000106,
          0.000124, 0.0001, 6.5e-5, 1.2e-5, 1.75e-5, 2.2e-5, 1.7e-5]

#--------------- U0 = 50 meV
# chemical potential
mu5 = [missing, missing, missing, -0.031, -0.030989, missing, -0.030941, missing, -0.030909, missing,
       -0.030808, missing, -0.030757, missing, -0.03069, missing, -0.030598, missing, -0.030588, missing,
       missing, -0.030194, missing, -0.029987, missing, missing, -0.029655, missing, -0.029295, -0.029079,
       -0.02891, -0.0287, -0.027963, -0.02695, -0.026367 , -0.025772, -0.024429]
# NN hopping
tnn5 = [missing, missing, missing, 0.00064, 0.000638, missing, 0.000632, missing, 0.000625, missing,
        0.000615, missing, 0.000602, missing, 0.000592, missing, 0.000582, missing, 0.000567, missing,
        missing, 0.000537, missing, 0.000497, missing, missing, 0.00046, missing, 0.000424, 0.000405,
        0.00039, 0.000373, 0.000315, 0.00026, 0.000233, 0.000206, 0.000157]
# NNN hopping
tnnn5 = [missing, missing, missing, 1.4e-5, 2.0e-5, missing, 2.1e-5, missing, 2.0e-5, missing,
         1.7e-5, missing, 1.7e-5, missing, 1.4e-5, missing, 1.65e-5, missing, 1.7e-5, missing,
         missing, 2.2e-5, missing, 1.3e-5, missing, missing, 1.1e-5, missing, 3.0e-6, 9.0e-6,
         9.0e-6, 9.0e-6, 6.0e-6, 2.0e-6, 5.0e-6, 3.0e-6, 8.0e-6]
# d-NNN hopping
tdnnn5 = [missing, missing, missing, 1.5e-5, 9.0e-6, missing, 1.4e-5, missing, 0.000152, missing,
          2.6e-5, missing, 5.8e-5, missing, 0.00019, missing, 4.1e-5, missing, 7.0e-5, missing,
          missing, 2.5e-5, missing, 1.2e-5, missing, missing, 3.0e-5, missing, 0.000176, 1.7e-5,
          2.0e-5, 2.0e-5, 1.3e-5, 2.0e-6, 2.0e-6, 6.0e-6, 2.0e-6] 

#--------------- U0 = 80 meV
# chemical potential
mu8 = [missing, missing, missing, -0.0703, -0.07014, missing, -0.07011, missing, missing, missing,
       -0.070028, missing, -0.069975, missing, -0.069925, missing, -0.06986, missing, missing, missing,
       missing, -0.069422, missing, missing, missing, missing, -0.069156, missing, -0.068867, missing,
       missing, missing, -0.06779, missing, -0.066672, missing, -0.065098]
# NN hopping
tnn8 = [missing, missing, missing, 0.000241, 0.000248, missing, 0.000246, missing, missing, missing,
        0.000241, missing, 0.000236, missing, 0.000235, missing, 0.000234, missing, missing, missing,
        missing, 0.000201, missing, missing, missing, missing, 0.000194, missing, 0.000185, missing,
        missing, missing, 0.000164, missing, 0.000113, missing, 7.9e-5]
# NNN hopping
tnnn8 = [missing, missing, missing, 4.0e-6, 5.0e-6, missing, 2.0e-6, missing, missing, missing,
         2.0e-6, missing, 8.0e-6, missing, 2.0e-6, missing, 5.0e-6, missing, missing, missing,
         missing, 5.2e-5, missing, missing, missing, missing, 4.0e-6, missing, 7.0e-6, missing,
         missing, missing, 1.7e-5, missing, 1.6e-5, missing, 2.3e-5] 
# d-NNN hopping
tdnnn8 = [missing, missing, missing, 1.7e-5, 1.0e-6, missing, 2.0e-6, missing, missing, missing,
          0.0, missing, 1.2e-5, missing, 9.0e-5, missing, 6.0e-6, missing, missing, missing,
          missing, 1.4e-5, missing, missing, missing, missing, 1.3e-5, missing, 2.0e-5, missing,
          missing, missing, 6.0e-6, missing, 1.0e-6, missing, 2.7e-5] 



## ---------------------- plotting 

pqplotlist = [(1, 8), (1, 3), (1, 2), (2, 3),(5, 6), (1, 1), (7,6), (4, 3),  (3, 2), (7, 4), (2, 1), (7, 3), (5, 2), (8, 3), (3, 1)]
phiplotlist = [pj/qj for (pj,qj) in pqplotlist]
labellist = ["$pj/$qj" for (pj,qj) in pqplotlist]

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