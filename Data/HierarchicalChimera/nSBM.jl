using Kuramodel
using Distributions
using Random
using JLD2
using LinearAlgebra

folderpath = "/mnt/lustre/scratch/inf/ec627/data/HierarchicalChimera/paper_data/raw_data/"

# job input
test = false
seedval = parse(Int, ARGS[1])
Random.seed!(seedval)

H_range = 0.00:0.01:1
for H in H_range
    # hSBM parameters
    n = [16, 8, 2]
    B = [16, 2, 1]
    N = prod(n)
    k = N / 3
    K = 1

    # build network model
    Random.seed!(seedval); A, P, p = SBMvar(n, B, H, k; K = K)

    # simulation parameters
    β = 0.05
    ω = repeat([1], N)
    if test
        sim_time = 1
    else
        sim_time = 35
    end
    Δt = 1e-4
    save_ratio = 1

    # simulation
    # lag parameter
    α = π/2 - β
        
    # construct lag matrix
    α_mat = zeros(N, N) # phase lag matrix
    for i in 1:N
        for j in i+1:N
            # add lag if different partition at layer 1
            if A[i,j] != 0 && P[1, i] != P[1, j]
                α_mat[i, j] = α
                α_mat[j, i] = α
            end
        end
    end

    # oscillators parameters
    noise_scale = 0

    # simulation parameters
    steps = (0.0+Δt):Δt:sim_time
    no_steps = length(steps)

    # storing parameters
    no_saves = round(Integer, no_steps / save_ratio)
    Random.seed!(seedval); global θ_now = rand(Uniform(-π, π), N)  # random init conditions
    θs = zeros(no_saves, N)
    θs[1, :] = θ_now

    global save_counter = 1

    for t in 2:no_steps
        
        # update phases
        θj_θi_mat = (repeat(θ_now',N) - repeat(θ_now',N)') - α_mat
        setindex!.(Ref(θj_θi_mat), 0.0, 1:N, 1:N) # set diagonal elements to zero 

        k1 = map(sum, eachrow(A .* sin.(θj_θi_mat)))
        global θ_now += Δt .* (ω + k1) + noise_scale*(rand(Normal(0,1),N))*sqrt(Δt)
        global save_counter += 1

        # save θ
        if save_counter % save_ratio == 0
            θs[round(Integer, save_counter / save_ratio), :] = θ_now
        end
        
    end


    # layer 1 data (modules)
    modules_KOP = zeros(B[1])
    modules_KOP_std = zeros(B[1])
    modules_macro_op = zeros(length(θs[:, 1]), B[1])

    for i in 1:B[1]
        modules_macro_op[:, i] = macro_op(θs[:, n[1]*(i-1)+1:n[1]*i])
        modules_KOP[i] = mean(modules_macro_op[:, i])
        modules_KOP_std[i] = std(modules_macro_op[:, i])
    end

    # layer 2 data (populations)
    if test
        relax = round(Integer, (1 / Δt) * 0.5)
    else
        relax = round(Integer, (1 / Δt) * 5)
    end
    mean_pop1 = round(mean(macro_op(θs[relax+1:end, 1:Integer(N/2)])), digits = 4)
    std_pop1 = round(std(macro_op(θs[relax+1:end, 1:Integer(N/2)])), digits = 4)
    mean_pop2 = round(mean(macro_op(θs[relax+1:end, Integer(N/2):N])), digits = 4)
    std_pop2 = round(std(macro_op(θs[relax+1:end, Integer(N/2):N])), digits = 4)

    params = Dict(
        "n" => n,
        "B" => B,
        "N" => N,
        "α" => α,
        "β" => β,
        "K" => K,
        "B" => B,
        "A" => A,
        "seedval" => seedval,
        "Δt" => Δt,
        "sim_time" => sim_time,
        "H" => H,
        "P" => P,
        "k" => k
    )

    # save results
    results = [θs, params]

    fileseed = "seed_" * string(seedval)
    filebeta = "_beta_" * string(β)
    fileH = "_H_" * string(H)
    filek = "_k_" * string(k)[1:4]

    filename = folderpath * fileseed * filebeta * fileH * filek

    if test
        println("end test: ", filename * ".jld2")
        save_object(filename * "_TEST.jld2", results)
        println("File saved correctly: $(filename)")
    else
        save_object(filename * ".jld2", results)
        println("File saved correctly: $(filename)")
    end
end