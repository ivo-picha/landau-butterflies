# job file that submits jobs of main_S.jl for sets of parameters

list_U0 = round.(collect(range(0.0010, 0.016, step = 0.0001)); digits=4)
#list_U0 = [0.0019, 0.004]
list_a = [5]
list_q = [50]

ip = collect(Base.product(list_U0, list_a, list_q))
param_list_tuple = reshape(ip, :)
param_list_str = replace.(string.(param_list_tuple), "(" => "[", ")" => "]")

n_cpus = 2 #number of cpus per job

folder_path = "/users/ivoga/lh/jobs"
output_msgs_path = "/users/ivoga/lh/msgs"


for (j,params) in enumerate(param_list_tuple)
    #check that an output hasn't already been generated for these parameters?? to be added
    println("submitting q_and_Ucrit.jl with parameters $params ($j of $(length(param_list_tuple)))")
    params_str = replace(string(params), "(" => "", ")" => "", ", " => "_", "," => "_")
    jobname = string("qS_", params_str, ".job")
    path_job = joinpath(folder_path, jobname)

    # create job file
    open(path_job, "w") do job
        #settings
        write(job, "#\$ -q rademaker \n")
        write(job, "#\$ -o $output_msgs_path -e $output_msgs_path \n")
        write(job, "#\$ -pe openmp $n_cpus \n")
        write(job, "#\$ -v OMP_NUM_THREADS=$n_cpus \n")
        write(job, "#\$ -v OMP_DYNAMIC=FALSE \n")
        write(job, "#\$ -v JULIA_NUM_THREADS=$n_cpus \n\n")

        #run file -------------------------------------------------------------------------------------------------------- OPTIONS GO BELOW ----------
        write(job, "cd /users/ivoga/lh/code/ \n") # change to project directory with .toml files
        write(job, "julia --project=. ./short/q_and_Ucrit.jl $(params[1]) $(params[2]) $(params[3]) \n")
    end

    run(`qsub $path_job`)   
end
