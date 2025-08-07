# job file that runs the main_D.jl for sets of parameters
# arguments should be in the format [p, q, U0_in_eV, a_in_Å, N_LLs, np-fillingfactor]
using Base.Threads

list_p = [1,2,3,4,5,6,7,8]
list_q = list_p
list_U0 = [0.05]
list_a = [50]
list_NLL = [20]
list_np = [1]
list_T = [1]

ip = collect(Base.product(list_p, list_q, list_U0, list_a, list_NLL, list_np, list_T))
param_list_tuple = reshape(ip, :)
param_list_str = replace.(string.(param_list_tuple), "(" => "[", ")" => "]")

@threads for j in eachindex(param_list_str)
    pj = param_list_tuple[j][1]
    qj = param_list_tuple[j][2]
    opt = param_list_str[j]
    if pj/qj == 1 && pj != 1
        continue
    end
    if pj/qj > 3 || pj/qj < 0.3
        continue
    end
    run(`julia ../mains/main_D.jl "$opt"`)
end

