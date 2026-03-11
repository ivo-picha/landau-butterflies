# job file that submits jobs of main_S.jl for sets of parameters

#list_U0 = round.(collect(range(0.003, 0.03, 8)); digits=4)
list_U0 = [0.05, 0.03]
list_p = collect(1:10)
list_q = list_p

list_used_phis = Tuple{Float32,Float32}[];

ip = collect(Base.product(list_U0, list_p, list_q))
param_list_tuple = reshape(ip, :)
param_list_str = replace.(string.(param_list_tuple), "(" => "[", ")" => "]")

n_cpus = 8 #number of cpus per job

folder_path = "/users/ivoga/lh/jobs"
output_msgs_path = "/users/ivoga/lh/msgs"


for (j,params) in enumerate(param_list_tuple)
    phi = round(Float32(params[2]/params[3]), digits=5)
    U0 = Float32(params[1])
    if (phi, U0) in list_used_phis
        continue
    elseif phi > 2.0 || phi < 0.2
        continue
    else
        push!(list_used_phis, (phi, U0))
    end
    LLmax = Int(round(30*25*params[1]/phi))
    #check that an output hasn't already been generated for these parameters?? to be added
    println("submitting main_t.jl with parameters $params and LLmax = $LLmax")
    params_str = replace(string(params), "(" => "", ")" => "", ", " => "_", "," => "_")
    jobname = string("t_", params_str, ".job")
    path_job = joinpath(folder_path, jobname)

    #create job file
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
        write(job, "julia --project=. ./mains/main_t.jl $(params[2]) $(params[3]) $(params[1]) $(LLmax)\n")
    end

    run(`qsub $path_job`)   
end
