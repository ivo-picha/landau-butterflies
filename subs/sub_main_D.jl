# job file that runs the main_D.jl for sets of parameters
# arguments should be in the format [p, q, U0_in_eV, a_in_Å, N_LLs, np-fillingfactor]

list_p = [1,2,3,4,5,6,7,8,9,10]
list_q = list_p
list_U0 = [0.015, 0,005]
list_a = [50]
list_NLL = [6,8,10,12,14,16]
list_np = [1]
list_T = [1]

ip = collect(Base.product(list_p, list_q, list_U0, list_a, list_NLL, list_np, list_T))
param_list_tuple = reshape(ip, :)
param_list_str = replace.(string.(param_list_tuple), "(" => "[", ")" => "]")

n_cpus = 1 #number of cpus per job

folder_path = "/users/ivoga/lh/jobs"
output_msgs_path = "/users/ivoga/lh/msgs"

for (j,params) in enumerate(param_list_str)

    #only create job in region that i want
    pj = param_list_tuple[j][1]
    qj = param_list_tuple[j][2]
    if pj/qj < 0.25 || pj/qj > 2
        continue
    elseif (pj/qj==0.5 && pj!=1) || (pj/qj==1 && pj !=1) || (pj/qj==1/3 && pj!=1)
        continue
    end

    #check that an output hasn't already been generated for these parameters?? to be added
    println("working on the set of parameters $params")
    params_str = replace(params, "[" => "", "]" => "", ", " => "_", "," => "_")
    jobname = string("D_", params_str, ".job")
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
        write(job, "julia -t $n_cpus ../mains/main_D.jl \"$params\" \n")
    end

    run(`qsub $path_job`)
end

