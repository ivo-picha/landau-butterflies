using Plots
using ImageMagick
using FileIO
using ColorTypes
using ImageCore

in_folder = "out_loc/wannier_states"
out_folder = "out_loc/gifs"

file_list = readdir(in_folder; join=true)
#sorted_file_list = sort(file_list, by = f -> parse(Float64, split(basename(f), "_")[3]))
files_quoted = join(["'$(f)'" for f in file_list], " ")

#frames = [load(file) for file in sorted_file_list]
#frames_wrapped = [ImageCore.Image(frame) for frame in frames]
#save(joinpath(out_folder, "spectrum_evolution_2BFs.gif"), frames; fps=6, loop=true)

cmd = "convert -delay 80 -loop 0 $files_quoted '$out_folder/wannier_states_circling.gif'"
println("Running: $cmd")
run(`sh -c $cmd`)