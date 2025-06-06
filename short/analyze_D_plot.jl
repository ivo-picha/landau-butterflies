include(joinpath(dirname(@__DIR__), "mods/densities.jl"))
using .Dens
include(joinpath(dirname(@__DIR__), "mods/params.jl"))
using .Params


using Plots
using NPZ
using ProgressMeter
using LaTeXStrings

folder_path = "/home/ivoga/Documents/PhD/Landau_Hofstadter/jl/data/mafalda/densities"

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

    return (sigmaR, dens_ratio, rot_score, parse(Int,p), parse(Int,q), parse(Float64,U0), parse(Int, NLL))
end

################ ADD HERE MANUALLY SOME PARAMETERS
# available NLLs
Nlist = [6,8,10,12, 14, 16]
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
        sR, dr, rs, p, q, U0, NLL = extract_data(fp)
        if U0 != 0.005
            continue
        end
        if (p/q==0.5 && p!=1) || (p/q==1 && p !=1) || (p/q==1/3 && p!=1) || (p/q==1/4 && p!=1) || (p/q==2/5 && p!=1)
            continue
        end
        x_axis_variable = p/q
        scatcolor = findfirst(x-> x==NLL, Nlist)
        scatter!(plt_sR, (x_axis_variable, sR), color = scatcolor, label = "")
        scatter!(plt_dr, (x_axis_variable, dr), color = scatcolor, label = "")
        scatter!(plt_rs, (x_axis_variable, rs), color = scatcolor, label = "")
    end
end


#create a legend with NLL
legend_plot = plot(legend=:outerright, grid=false, framestyle=:none)
for j in eachindex(Nlist)
    scatter!(legend_plot, (NaN, NaN), color=j, label="$(Nlist[j]+1) LLs")
end


name_plots_base = "varyphi_0to2_U0.005_6Ns_"
savefig(plt_sR, joinpath(plots_save_path, string(name_plots_base, "sR.png")))
savefig(plt_dr, joinpath(plots_save_path, string(name_plots_base, "dr.png")))
savefig(plt_rs, joinpath(plots_save_path, string(name_plots_base, "rs.png")))
savefig(legend_plot, joinpath(plots_save_path, string(name_plots_base, "leg.png")))
#plt_sR
#plot(plt_sR, plt_dr, plt_rs; layout=(3,1), size=(400,1200), link = :x)