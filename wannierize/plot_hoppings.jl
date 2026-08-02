using Plots
using Measures

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
# Each tuple contains:
# (
#   p,
#   q,
#   index_to_m,
#   mixed,
#   mean_abs_t_dNNN,
#   residual_phase_over_pi,
#   coherence
# )
dNNN_results_20meV = [
    (0, 1, [0],                    false, 1.4311e-05, +0.0000, 1.000),
    (1, 7, [6, 5, 1, 0, 4, 2, 3], true,  5.3041e-05, +0.0227, 0.986),
    (1, 6, [1, 5, 4, 0, 2, 3],    true,  1.5747e-04, -0.1017, 0.739),
    (1, 5, [0, 4, 3, 1, 2],       true,  9.3446e-05, -0.0314, 0.983),
    (1, 4, [0, 3, 1, 2],          false, 1.3184e-04, -0.0349, 0.984),
    (2, 7, [3, 6, 0, 2, 4, 5, 1], false, 4.1450e-04, +0.0096, 0.725),
    (1, 3, [0, 2, 1],             false, 1.9721e-04, -0.0387, 0.987),
    (2, 5, [2, 4, 0, 1, 3],       true,  3.5530e-04, -0.1104, 0.622),
    (3, 7, [4, 2, 0, 6, 5, 3, 1], true,  4.7034e-04, -0.0774, 0.925),
    (1, 2, [1, 0],                false, 3.1739e-04, +0.0527, 0.990),
    (4, 7, [3, 5, 1, 6, 0, 2, 4], false, 4.3782e-04, -0.0698, 0.961),
    (3, 5, [1, 3, 0, 4, 2],       false, 4.9499e-04, -0.0993, 0.969),
    (2, 3, [1, 0, 2],             false, 6.0618e-04, -0.0951, 0.907),
    (5, 7, [1, 5, 4, 2, 0, 6, 3], false, 4.0292e-04, -0.0830, 0.981),
    (3, 4, [2, 1, 3, 0],          false, 4.4687e-04, -0.0832, 0.804),
    (4, 5, [2, 3, 1, 4, 0],       true,  3.6578e-04, -0.0494, 0.862),
    (5, 6, [2, 3, 0, 4, 5, 1],    true,  3.2541e-04, -0.0303, 0.861),
    (6, 7, [3, 4, 2, 0, 5, 1, 6], true,  3.0965e-04, -0.0446, 0.989),
    (1, 1, [0],                    false, 2.1100e-04, +0.0029, 1.000),
    (7, 6, [2, 3, 4, 1, 5, 0],    true,  5.7515e-05, +0.1791, 0.970),
    (6, 5, [2, 1, 3, 0, 4],       true,  4.1788e-05, +0.2517, 0.905),
    (5, 4, [2, 1, 0, 3],          true,  6.6505e-05, +0.8915, 0.630),
    (4, 3, [1, 2, 0],             false, 7.8863e-05, +0.5896, 0.909),
    (7, 5, [3, 1, 0, 4, 2],       false, 5.9275e-05, +0.8378, 0.973),
    (3, 2, [0, 1],                false, 1.0672e-04, +0.5522, 0.965),
    (8, 5, [2, 4, 0, 3, 1],       false, 9.0709e-05, +0.9206, 0.988),
    (5, 3, [0, 2, 1],             false, 9.0217e-05, +0.9306, 1.000),
    (7, 4, [0, 3, 2, 1],          false, 8.8271e-05, +0.9573, 0.999),
    (2, 1, [0],                    false, 6.5000e-05, +1.0000, 1.000),
    (7, 3, [2, 0, 1],             false, 1.1053e-05, -0.4062, 0.640),
    (5, 2, [0, 1],                false, 1.5653e-05, -0.2611, 0.733),
    (8, 3, [1, 2, 0],             false, 1.6778e-05, -0.0848, 0.979),
    (3, 1, [0],                    false, 1.7000e-05, +0.0111, 1.000),
]
dnnn_phi_20 = Float64[]
dnnn_abs_20 = Float64[]
dnnn_arg_20 = Float64[]
for dnnn in dNNN_results_20meV
    phi = dnnn[1]/dnnn[2]
    yarg = dnnn[6]
    yabs = dnnn[5]
    if dnnn[end]>0.7 #|| !dnnn[4]
        push!(dnnn_phi_20,phi)
        push!(dnnn_abs_20,yabs)
        push!(dnnn_arg_20,yarg)
    end
end

#--------------- U0 = 50 meV
# chemical potential
mu5 = [-0.031030036,
       missing, missing, -0.030994, -0.031, -0.030989, missing, -0.030941, missing, -0.030909, missing,
       -0.030808, missing, -0.030757, missing, -0.03069, missing, -0.030598, missing, -0.030588, missing,
       missing, -0.030194, missing, -0.029987, missing, missing, -0.029655, missing, -0.029295, -0.029079,
       -0.02891, -0.0287, -0.027963, -0.02695, -0.026367 , -0.025772, -0.024429]
# NN hopping
tnn5 = [0.0006277,
        missing, missing, 0.000645, 0.00064, 0.000638, missing, 0.000632, missing, 0.000625, missing,
        0.000615, missing, 0.000602, missing, 0.000592, missing, 0.000582, missing, 0.000567, missing,
        missing, 0.000537, missing, 0.000497, missing, missing, 0.00046, missing, 0.000424, 0.000405,
        0.00039, 0.000373, 0.000315, 0.00026, 0.000233, 0.000206, 0.000157]
# NNN hopping
tnnn5 = [1.5045e-7,
         missing, missing, 2.0e-5, 1.4e-5, 2.0e-5, missing, 2.1e-5, missing, 2.0e-5, missing,
         1.7e-5, missing, 1.7e-5, missing, 1.4e-5, missing, 1.65e-5, missing, 1.7e-5, missing,
         missing, 2.2e-5, missing, 1.3e-5, missing, missing, 1.1e-5, missing, 3.0e-6, 9.0e-6,
         9.0e-6, 9.0e-6, 6.0e-6, 2.0e-6, 5.0e-6, 3.0e-6, 8.0e-6]
# d-NNN hopping
dNNN_results_50meV = [
    (0, 1, [0],                   false, 1.474e-6,  0.000, 1.00),
    (1, 7, [6, 1, 5, 0, 4, 2, 3], true,  4.125e-6,  0.037, 0.84),
    (1, 6, [0, 5, 1, 4, 2, 3],    true,  4.875e-6, -0.025, 0.99),
    (1, 5, [0, 4, 3, 1, 2],       true,  7.327e-6, -0.020, 0.88),
    (1, 4, [0, 3, 1, 2],          false, 8.882e-6, -0.027, 0.99),
    (2, 7, [3, 6, 0, 2, 4, 1, 5], false, 1.032e-4,  0.008, 0.52),
    (1, 3, [0, 2, 1],             false, 1.432e-5, -0.029, 0.99),
    (2, 5, [2, 0, 4, 1, 3],       true,  1.004e-4, -0.286, 0.38),
    (3, 7, [2, 4, 0, 6, 5, 1, 3], true,  6.013e-5, -0.246, 0.63),
    (1, 2, [1, 0],                false, 2.584e-5,  0.031, 1.00),
    (4, 7, [3, 5, 1, 6, 0, 4, 2], false, 6.500e-5, -0.082, 0.91),
    (3, 5, [1, 3, 4, 0, 2],       false, 8.849e-5, -0.088, 0.80),
    (2, 3, [1, 0, 2],             true,  7.798e-5, -0.150, 0.58),
    (5, 7, [5, 1, 4, 2, 0, 3, 6], true,  1.834e-5, -0.279, 0.73),
    (3, 4, [1, 2, 0, 3],          false, 9.727e-5, -0.207, 0.49),
    (5, 6, [2, 3, 0, 4, 5, 1],    true,  4.262e-5,  0.083, 0.51),
    (6, 7, [3, 2, 4, 0, 1, 5, 6], true,  4.356e-5, -0.054, 0.61),
    (1, 1, [0],                   false, 2.500e-5,  0.000, 1.00),
    (7, 6, [3, 2, 4, 0, 1, 5],    true,  1.507e-5,  0.251, 0.85),
    (4, 3, [1, 2, 0],             false, 2.563e-5,  0.512, 1.00),
    (7, 5, [1, 3, 0, 2, 4],       true,  3.814e-5,  0.846, 0.45),
    (3, 2, [1, 0],                false, 3.877e-5,  0.495+0.2, 1.00),
    (8, 5, [2, 4, 0, 1, 3],       false, 1.465e-5,  0.907, 0.99),
    (5, 3, [0, 2, 1],             false, 1.548e-5,  0.926, 0.96),
    (7, 4, [3, 0, 2, 1],          false, 1.557e-5,  0.967, 0.98),
    (2, 1, [0],                   false, 1.300e-5,  1.000, 1.00),
    (7, 3, [0, 2, 1],             false, 1.491e-6,  0.082, 0.13),
    (5, 2, [0, 1],                false, 1.949e-6, -0.150, 0.55),
    (8, 3, [1, 0, 2],             false, 3.659e-6, -0.076, 0.98),
    (3, 1, [0],                   false, 2.000e-6, -0.829+1, 1.00),
]
dnnn_phi_50 = Float64[]
dnnn_abs_50 = Float64[]
dnnn_arg_50 = Float64[]
for dnnn in dNNN_results_50meV
    phi = dnnn[1]/dnnn[2]
    yarg = dnnn[6]
    yabs = dnnn[5]
    if dnnn[end]>0.7 #|| !dnnn[4]
        push!(dnnn_phi_50,phi)
        push!(dnnn_abs_50,yabs)
        push!(dnnn_arg_50,yarg)
    end
end


#--------------- U0 = 80 meV
# chemical potential
mu8 = [-0.0701928,
       missing, missing, -0.07004, -0.0703, -0.07014, missing, -0.07011, missing, -0.070094, missing,
       -0.070028, missing, -0.069975, missing, -0.069925, missing, -0.06986, missing, missing, missing,
       missing, -0.069422, missing, missing, -0.06936, missing, -0.069156, missing, -0.068867, missing,
       missing, -0.068451, -0.06779, -0.067127, -0.066672, -0.066205, -0.065098]
# NN hopping
tnn8 = [0.0002035504457158282,
        missing, missing, missing, 0.000241, 0.000248, missing, 0.000246, missing, 0.000244, missing,
        0.000241, missing, 0.000236, missing, 0.000235, missing, 0.000234, missing, missing, missing,
        missing, 0.000201, missing, missing, missing, missing, 0.000194, missing, 0.000185, missing,
        missing, 0.000166, 0.000164, 0.000125, 0.000113, 0.0001, 7.9e-5]
# NNN hopping
tnnn8 = [4.32333950e-5,
         missing, missing, missing, 4.0e-6, 5.0e-6, missing, 2.0e-6, missing, 2.0e-6, missing,
         2.0e-6, missing, 8.0e-6, missing, 2.0e-6, missing, 5.0e-6, missing, missing, missing,
         missing, 5.2e-5, missing, missing, missing, missing, 4.0e-6, missing, 7.0e-6, missing,
         missing, 1.0e-6, 1.7e-5, 1.0e-6, 1.6e-5, 5.0e-6, 2.3e-5] 
# d-NNN hopping
dNNN_results_80meV = [
    (0, 1, [0],                    false, 2.8144e-07, +0.0000, 1.000),
    (1, 7, [5, 4, 2, 6, 1, 3, 0], true,  1.2004e-05, -0.0216, 0.455),
    (1, 6, [0, 4, 1, 3, 5, 2],    true,  1.9266e-05, +0.2445, 0.850),
    (1, 5, [1, 0, 4, 2, 3],       true,  2.3752e-06, -0.3035, 0.953),
    (1, 4, [0, 3, 1, 2],          false, 9.9431e-07, -0.0186, 0.995),
    (1, 3, [0, 2, 1],             false, 1.9990e-06, -0.0214, 0.995),
    (2, 5, [2, 4, 0, 1, 3],       false, 3.8455e-05, +0.0083, 0.489),
    (1, 2, [1, 0],                false, 2.9899e-06, +0.0208, 0.998),
    (4, 7, [3, 5, 1, 6, 0, 4, 2], false, 1.5043e-05, -0.0427, 0.860),
    (3, 5, [1, 3, 0, 4, 2],       false, 2.2686e-05, -0.0525, 0.715),
    (2, 3, [1, 0, 2],             false, 4.7757e-05, -0.1930, 0.620),
    (5, 7, [1, 5, 4, 2, 6, 0, 3], false, 1.4533e-05, -0.1314, 0.456),
    (3, 4, [2, 1, 0, 3],          false, 3.0432e-05, -0.3167, 0.320),
    (4, 5, [2, 4, 1, 3, 0],       true,  8.4261e-05, -0.4679, 0.875),
    (1, 1, [0],                    false, 1.4000e-05, +0.9637-1, 1.000),
    (8, 7, [3, 5, 4, 0, 6, 2, 1], true,  1.2384e-04, -0.8272+1, 0.843),
    (7, 6, [3, 2, 4, 0, 1, 5],    true,  2.9066e-06, +0.1092, 0.786),
    (6, 5, [2, 1, 3, 4, 0],       true,  1.1620e-05, +0.5469, 0.938),
    (5, 4, [2, 1, 3, 0],          false, 7.6069e-06, +0.4579, 0.963),
    (4, 3, [1, 2, 0],             false, 1.0254e-05, +0.4810, 0.992),
    (7, 5, [1, 3, 4, 0, 2],       false, 3.6680e-06, +0.5371, 0.897),
    (3, 2, [0, 1],                false, 1.4493e-05, +0.4739, 0.998),
    (8, 5, [2, 0, 4, 1, 3],       false, 2.5917e-06, +0.9004, 0.953),
    (5, 3, [0, 1, 2],             true,  1.8888e-05, -0.3162+1, 0.941),
    (7, 4, [0, 3, 2, 1],          false, 3.6600e-06, -0.9539+2, 0.872),
    (2, 1, [0],                    false, 6.0000e-06, -0.6863+1.6, 1.000),
    (7, 3, [0, 2, 1],             false, 9.9946e-07, +0.4440, 0.768),
    (5, 2, [1, 0],                false, 9.9979e-07, +0.5562, 0.978),
    (8, 3, [1, 0, 2],             false, 8.2848e-07, -0.1909, 1.000),
    (3, 1, [0],                    true,  2.7000e-05, -0.5586, 1.000),
]
dnnn_phi_80 = Float64[]
dnnn_abs_80 = Float64[]
dnnn_arg_80 = Float64[]
for dnnn in dNNN_results_80meV
    phi = dnnn[1]/dnnn[2]
    yarg = dnnn[6]
    yabs = dnnn[5]
    if dnnn[end]>0.7 && !(dnnn[4])
        push!(dnnn_phi_80,phi)
        push!(dnnn_abs_80,yabs)
        push!(dnnn_arg_80,yarg)
    end
end



## y NN hopping phases, multiples of pi; for U = 50 meV only
ph_list = [
    # q = 1
    (1, 1, [0.012694]),
    (2, 1, [0.990117]),
    (3, 1, [0.002523]),

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
    dpi = 120,
    ms = 8,
)
for n in 1:5
    r = 0.01 + n / 4
    plot!(
        circles,
        r .* cos.(θ),
        r .* sin.(θ),
        label="",
        color = :black,
        lw=2
    )
end

function plot_phases!(plot::Plots.Plot, set::Tuple{Int64,Int64,Vector})
        p, q, pts = set
        r = 0.01 + q / 4
        ptsxy = [(r*cos(pt*π), r*sin(pt*π)) for pt in pts]
        scatter!(plot, ptsxy, color = (iseven(p) ? :red : :blue), label = "", ms=7)
end


gp1 = π.*[0,1]
gp2 = π.*[1/2,3/2]
gp3 = π.*[0,1/3,2/3,4/3,5/3,1]
gp4 = π.*[1/4,3/4,5/4,7/4]
gp5 = π.*[0,1/5,2/5,3/5,4/5,1,6/5,7/5,8/5,9/5]

function plot_line_phases!(plt::Plots.Plot, phase::Number, q::Integer)
        r1 = 0.01 + (q-1) / 4
        r2 = 0.01 + q / 4
        x1, y1 = r1*cos(phase), r1*sin(phase)
        x2, y2 = r2*cos(phase), r2*sin(phase)
        trange = range(0,1,20)
        data = [(x1+(x2-x1)*t, y1+(y2-y1)*t) for t in trange]
        plot!(plt, data, label = "", lw = 2, style = :dot, color = :gray)
end

for phase in gp5
        plot_line_phases!(circles, phase, 5)
end
for phase in gp4
        plot_line_phases!(circles, phase, 4)
end
for phase in gp3
        plot_line_phases!(circles, phase, 3)
end
for phase in gp2
        plot_line_phases!(circles, phase, 2)
end
for phase in gp1
        plot_line_phases!(circles, phase, 1)
end

for set in ph_list
        plot_phases!(circles, set)
end

display(circles)
savefig(circles,"/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/phases.pdf")




## ---------------------- plotting 

pqplotlist = [(1, 8), (1, 3), (1, 2), (3,4), (1, 1), (7,6), (4, 3),  (3, 2), (7, 4), (2, 1), (7, 3), (5, 2), (8, 3), (3, 1)]
phiplotlist = [pj/qj for (pj,qj) in pqplotlist]
labellist = ["$pj/$qj" for (pj,qj) in pqplotlist]
pushfirst!(phiplotlist, 0)
pushfirst!(labellist, "0")

philist = [pj/qj for (pj,qj) in pqlist]

plt_opt = (
    framestyle = :box,
    xaxis = "ϕ = p/q",
    xlims = (0,3),
    size = (400,320),
    xticks = (phiplotlist,labellist),
    markers = :diamond,
    ms = 6,
    legend = :topright,
    legendfont = 11,
    guidefont = 13
)

plot_mu = scatter(philist,mu2.*1000, label= "";plt_opt...)
scatter!(plot_mu,philist,mu5.*1000, label= "";plt_opt...)
scatter!(plot_mu,philist,mu8.*1000, label= "";plt_opt...)
yaxis!(plot_mu, "μ [meV]")

plot_tnn = scatter(philist,tnn2.*1000, label= "U₀ = 20 meV";plt_opt...)
scatter!(plot_tnn,philist,tnn5.*1000, label= "U₀ = 50 meV";plt_opt...)
scatter!(plot_tnn,philist,tnn8.*1000, label= "U₀ = 80 meV";plt_opt...)
yaxis!(plot_tnn, "|t| [meV]")

plot_tnnn = scatter(philist,tnnn2.*1000, label= "";plt_opt...)
scatter!(plot_tnnn,philist,tnnn5.*1000, label= "";plt_opt...)
scatter!(plot_tnnn,philist,tnnn8.*1000, label= "";plt_opt...)
yaxis!(plot_tnnn, "|t'| [meV]")

plot_tdnnn = scatter(dnnn_phi_20,dnnn_abs_20.*1000, label= "";plt_opt...)
scatter!(plot_tdnnn,dnnn_phi_50,dnnn_abs_50.*1000, label= "";plt_opt...)
scatter!(plot_tdnnn,dnnn_phi_80,dnnn_abs_80.*1000, label= "";plt_opt...)
yaxis!(plot_tdnnn, "|t''| [meV]")

plot_mu_overlay = scatter(philist,(mu2.-minimum(skipmissing(mu2))).*1000, label= "";plt_opt...)
scatter!(plot_mu_overlay,philist,(mu5.-minimum(skipmissing(mu5))).*1000, label= "";plt_opt...)
scatter!(plot_mu_overlay,philist,(mu8.-minimum(skipmissing(mu8))).*1000, label= "";plt_opt...)
yaxis!(plot_mu_overlay, "μ - μ₀ [meV]")

plot_tnn_log = scatter(philist,log.(tnn2.*1000), label= "";plt_opt...)
scatter!(plot_tnn_log,philist,log.(tnn5.*1000), label= "";plt_opt...)
scatter!(plot_tnn_log,philist,log.(tnn8.*1000), label= "";plt_opt...)
yaxis!(plot_tnn_log, "log |t| [meV]")

plot_tdnnn_arg = scatter(dnnn_phi_20,dnnn_arg_20, label= "";plt_opt...)
scatter!(plot_tdnnn_arg,dnnn_phi_50,dnnn_arg_50, label= "";plt_opt...)
scatter!(plot_tdnnn_arg,dnnn_phi_80,dnnn_arg_80, label= "";plt_opt...)
yticks!(plot_tdnnn_arg, [0,0.5,1])
yaxis!(plot_tdnnn_arg, "arg(t'')/π")

plt_combi = plot(plot_mu, plot_tnn, plot_tnnn, plot_tdnnn, plot_tdnnn_arg, plot_mu_overlay, layout = (3,2), link=:x,
        size = (1.6*650,690), margins = 3mm, dpi = 250
)

savefig(plt_combi,"/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/hops_all2.pdf")





## extra plots, universal?
plot_mu_overlay = scatter(philist,(mu2.-minimum(skipmissing(mu2))).*1000, label= "";plt_opt...)
scatter!(plot_mu_overlay,philist,(mu5.-minimum(skipmissing(mu5))).*1000, label= "";plt_opt...)
scatter!(plot_mu_overlay,philist,(mu8.-minimum(skipmissing(mu8))).*1000, label= "";plt_opt...)
yaxis!(plot_mu_overlay, "μ - μ₀ [meV]")

plot_tnn_log = scatter(philist,log.(tnn2.*1000), label= "";plt_opt...)
scatter!(plot_tnn_log,philist,log.(tnn5.*1000), label= "";plt_opt...)
scatter!(plot_tnn_log,philist,log.(tnn8.*1000), label= "";plt_opt...)
yaxis!(plot_tnn_log, "log |t| [meV]")

plot_combi_extra = plt_combi = plot(plot_mu_overlay, plot_tnn_log, layout = (2,1),
        size = (1.6*300,500), margins = 6mm, dpi = 120
)


savefig(plot_combi_extra,"/home/ivoga/Documents/PhD/Landau_Hofstadter/jl2/out_loc/hops_abs_extra.pdf")