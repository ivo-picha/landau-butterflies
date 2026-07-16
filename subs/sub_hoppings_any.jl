
list_U0 = [0.02, 0.05, 0.08]
list_a = [5]
# list_LLmax = [625]
list_p = [1,2,3,4,5,6,7,8]
list_q = copy(list_p)
list_NXY = [64]
list_Nuc = [40]

ip = collect(Base.product(list_U0, list_a, list_p, list_q, list_NXY, list_Nuc))
param_list_tuple = reshape(ip, :)
param_list_str = replace.(string.(param_list_tuple), "(" => "[", ")" => "]")

n_cpus = 4 #number of cpus per job

folder_path = "/users/ivoga/lh/jobs"
output_msgs_path = "/users/ivoga/lh/msgs"

past_phi_U = [(0.0,0.0)]

for (j,params) in enumerate(param_list_tuple)

    phi = Float64(params[3]/params[4])
    U0 = Float64(params[1])

    if in((phi,U0),past_phi_U)
        continue
    elseif phi > 3
        continue
    end

    LLmax = clamp( Int64(round((25/phi) * (U0 * 20))), 10, 120)

    println(LLmax)

    #check that an output hasn't already been generated for these parameters?? to be added
    println("submitting hoppings_any.jl with parameters $params ($j of $(length(param_list_tuple)))")
    params_str = replace(string(params), "(" => "", ")" => "", ", " => "_", "," => "_")
    jobname = string("hops_", params_str, ".job")
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
        write(job, "julia --project=. ./wannierize/hoppings_any.jl $(params[3]) $(params[4]) $(params[2]) $(params[1]) $LLmax $(params[5]) $(params[6]) \n")
    end

    run(`qsub $path_job`)   

    push!(past_phi_U, (phi,U0))
end
