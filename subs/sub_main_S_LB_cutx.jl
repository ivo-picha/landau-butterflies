# job file that runs the main_SCWC_fixE_vp_LB.jl for sets of parameters
# arguments should be in the format [start_phi_inv, end_phi_inv, U0_in_eV, a_in_Å, p_bands, N_LLs, gap_factor]
# cuts the x axis and breaks it down into different jobs

Nxjobs = 14;

startphi = 0.05
endphi = 1.0

U0 = 0.03
a = 50
q = 120
Nmin = 14
gf = 0.05

phis_bd = round.(collect(range(startphi,endphi,Nxjobs+1)); digits = 4)
# make the jobs non-overlapping at the boundaries
plist = round.(phis_bd .* q)
endphis_p = round.((plist.-1)./q ; digits = 4)

#assess how many LLs are needed for each job to stay at a fixed max energy more or less 
NmaxLL(phi_n) = Int(round(endphi*(Nmin+0.5)/phi_n - 0.5))
Nlist = NmaxLL.(phis_bd)

param_list_tuple = Tuple[];
for j in 1:Nxjobs
    push!(param_list_tuple, (phis_bd[j], endphis_p[j+1], U0, a, q, Nlist[j+1], gf))
end
param_list_str = replace.(string.(param_list_tuple), "(" => "[", ")" => "]")

n_cpus = 1 #number of cpus per job

folder_path = "/users/ivoga/lh/jobs"
output_msgs_path = "/users/ivoga/lh/msgs"

for params in param_list_str
    #check that an output hasn't already been generated for these parameters?? to be added
    println("working on the set of parameters $params")
    params_str = replace(params, "[" => "", "]" => "", ", " => "_", "," => "_")
    jobname = string("S_LB_cut_", params_str, ".job")
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
        write(job, "julia -t $n_cpus ../mains/main_SCWC_fixE_vp_LB.jl \"$params\" \n")
    end

    run(`qsub $path_job`)    # could be replaced by $jobname if working in same folder
end

