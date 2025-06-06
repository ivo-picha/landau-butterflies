include(joinpath(dirname(@__DIR__), "mods/densities.jl"))
using .Dens
include(joinpath(dirname(@__DIR__), "mods/params.jl"))
using .Params


using Plots
using NPZ
using ProgressMeter
using LaTeXStrings

folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/mafalda/dens6_compare_np/"

function extract_data(file_path::String)
    if splitext(file_path)[2] != ".npz"
        error("File is not in .npz format!")
    end

    data = npzread(file_path)
    data_x = data["x"]
    data_y = data["y"]
    data_z = data["z"]

    file_name = basename(file_path)
    pd = Params.extract_params(file_name)

    p = pd[:p]
    q = pd[:q]
    np = parse(Float64, pd[:np])
    U0 = pd[:U0]
    a_in_angstr = pd[:a_in_angstr]
    a = parse(Float64,a_in_angstr) * 1e-10
    NLL = pd[:NLL]
    TK = pd[:TK]

    dxn, dyn, dzn = Dens.square_interpolate_d(data_x, data_y, data_z)
    
    sigmaR = Dens.get_sigmaR(dxn./a, dyn./a, dzn)
    dens_ratio = Dens.get_dens_ratio(dzn)
    rot_score = Dens.rotational_symmetry_score(dzn)

    return (sigmaR, dens_ratio, rot_score, parse(Int,p), parse(Int,q), parse(Float64,U0), parse(Int, NLL), np)
end

################ ADD HERE MANUALLY SOME PARAMETERS
# available NLLs
nlist = [0.9, 1.0, 1.1, 1.5, 2.0]#collect(0.005:0.003:0.03)
xax = L"\phi"
xl = (0.0,2.0)
plots_save_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/plots/local/random"

plt_opts = (framestyle=:box, size = (400,400), xlims = xl) #aspect_ratio=1.08*(xl[2]-xl[1])
plt_sR = plot(yaxis=L"\sigma_R";plt_opts...)
plt_dr = plot(yaxis=L"ρ_{max}/ρ_{bond}";plt_opts...)
plt_rs = plot(yaxis="Rot. corr.";plt_opts...)

# add x axis label
xaxis!(plt_sR, xax)
xaxis!(plt_dr, xax)
xaxis!(plt_rs, xax)

@showprogress for file in readdir(folder_path)
    fp = joinpath(folder_path, file)
    if isfile(fp)
        sR, dr, rs, p, q, U0, NLL, np = extract_data(fp)
        # if (p/q==0.5 && p!=1) || (p/q==1/3 && p!=1) || (p/q==1/4 && p!=1) || (p/q==2/5 && p!=1) || (p<=8) 
        #     continue
        # end
        # if q!=25
        #     continue
        # end
        x_axis_variable = p/q
        scatcolor = findfirst(x-> x==np, nlist)
        scatter!(plt_sR, (x_axis_variable, sR), color = scatcolor, label = "", ms = 3, markerstrokewidth = 0)
        scatter!(plt_dr, (x_axis_variable, dr), color = scatcolor, label = "", ms = 3, markerstrokewidth = 0)
        scatter!(plt_rs, (x_axis_variable, rs), color = scatcolor, label = "", ms = 3, markerstrokewidth = 0)
    end
end


#create a legend with NLL
legend_plot = plot(legend=:outerright, grid=false, framestyle=:none, size = (400,400))
for j in eachindex(nlist)
    scatter!(legend_plot, (NaN, NaN), color=j, label="nₑ=$(nlist[j])", markerstrokewidth = 0)
end

# add vlines on gaps (manual)
# xlines = [0.4,0.6,0.8,0.9,1.1,1.2,1.8,1.9]
# xlines_style = (label = "", color = :red, linestyle=:dash, alpha = 0.7)
# vline!(plt_sR, xlines; xlines_style...)
# vline!(plt_rs, xlines; xlines_style...)
# vline!(plt_dr, xlines; xlines_style...)


name_plots_base = "varyphi_0.2to2.0_U0.01_multinp_"
savefig(plt_sR, joinpath(plots_save_path, string(name_plots_base, "sR.png")))
savefig(plt_dr, joinpath(plots_save_path, string(name_plots_base, "dr.png")))
savefig(plt_rs, joinpath(plots_save_path, string(name_plots_base, "rs.png")))
savefig(legend_plot, joinpath(plots_save_path, string(name_plots_base, "zleg.png")))




# scatter!(plt_sR, (1.0,0.3830420934026403), label = "", ms = 1.5)
# plt_sR
#plot(plt_sR, plt_dr, plt_rs; layout=(3,1), size=(400,1200), link = :x)