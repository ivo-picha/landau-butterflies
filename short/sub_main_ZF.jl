# job file that runs the main_Dsm.jl for sets of parameters
# arguments should be in the format [p, q, U0_in_eV, a_in_Å, N_LLs, np-fillingfactor]


list_U0 = [0.01, 0.1, 1, 5, 10, 100]
list_nu = [1, 2, 3, 4]

ip = collect(Base.product(list_U0, list_nu))
param_list_tuple = reshape(ip, :)
param_list_str = replace.(string.(param_list_tuple), "(" => "[", ")" => "]")

n_cpus = 2 #number of cpus per job

folder_path = "/users/ivoga/lh-short/jobs"
output_msgs_path = "/users/ivoga/lh-short/msgs"

for params in param_list_str
    #check that an output hasn't already been generated for these parameters?? to be added
    println("working on the set of parameters $params")
    params_str = replace(params, "[" => "", "]" => "", ", " => "_", "," => "_")
    jobname = string("ZF_", params_str, ".job")
    path_job = joinpath(folder_path, jobname)

    # create job file
    open(path_job, "w") do job
        #settings
        write(job, "#\$ -q rademaker \n")
        write(job, "#\$ -o $output_msgs_path -e $output_msgs_path \n")
        write(job, "#\$ -pe openmp $n_cpus \n")
        write(job, "#\$ -v OMP_NUM_THREADS=$n_cpus \n")
        write(job, "#\$ -v OMP_DYNAMIC=FALSE \n\n")

        #run file
        write(job, "julia main_ZF.jl \"$params\" \n")
    end

    run(`qsub $path_job`)
end