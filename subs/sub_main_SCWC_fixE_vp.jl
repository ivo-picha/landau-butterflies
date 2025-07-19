# job file that runs the main_SCWC.jl for sets of parameters
# arguments should be in the format [start_phi_inv, end_phi_inv, U0_in_eV, a_in_Å, p_bands, N_LLs, gap_factor]

list_startphi = [0.25]
list_endphi = [2.0]
# list_U0 = collect(range(0.0010, 0.0645, step = 0.0005))
list_U0 = [0.05, 0.01, 0.015, 0.02]
list_a = [50]
list_q = [360]
list_Nmin = [3]
list_gf = [0.08]

ip = collect(Base.product(list_startphi, list_endphi, list_U0, list_a, list_q, list_Nmin, list_gf))
param_list_tuple = reshape(ip, :)
param_list_str = replace.(string.(param_list_tuple), "(" => "[", ")" => "]")

n_cpus = 16 #number of cpus per job

folder_path = "/users/ivoga/lh/jobs"
output_msgs_path = "/users/ivoga/lh/msgs"

for params in param_list_str
    #check that an output hasn't already been generated for these parameters?? to be added
    println("working on the set of parameters $params")
    params_str = replace(params, "[" => "", "]" => "", ", " => "_", "," => "_")
    jobname = string("SCWC_fixE_vp_", params_str, ".job")
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

        #run file
        write(job, "julia -t $n_cpus ../mains/main_SCWC_fixE_vp.jl \"$params\" \n")
    end

    run(`qsub $path_job`)    # could be replaced by $jobname if working in same folder
end

