using Kuramodel
using JLD2
using Random
using Distributions
using LinearAlgebra
using MCIntegration
using Graphs

seedval = 123
Random.seed!(seedval)

test = false

# hyperbolic network parameters
N = 256
k = 8
if test
    γ_tests = [2.1, 2.9]
else
    γ_tests = collect(2.01:0.01:2.99)
end

for γ in γ_tests

    ζ = 1
    if test
        T = 0.1
    else
        T = parse(Int, ARGS[1]) / 100
    end
    graphtype = :softHRGG
    if graphtype == :HRGG	
        HRGGnetwork = hyperbolic_graph(:HRGG, N, k, γ, 0.0; ζ = ζ, ϵ = ϵ, max_iter = max_iter)
    elseif graphtype == :softHRGG
        HRGGnetwork = hyperbolic_graph(:softHRGG, N, k, γ, T; ζ = ζ)
    end
    A = HRGGnetwork[1]

    # oscillators parameters
    nat_f = 40  # Hz
    ω = repeat([2 * π * nat_f], N)
    σ = repeat([1], N)
    K = 8 # global coupling

    # simulation parameters
    Δt = 1e-4 # 0.1 ms
    A = A .* Δt .* K # rescaled by the time step and coupling already
    if test
        sim_time = 0.1  # seconds
    else
        sim_time = 20  # seconds
    end
    steps = (0.0 + Δt):Δt:sim_time
    no_steps = length(steps)

    # delay Matrix
    τ = zeros(N, N)
    τ_global = 4 * 1e-3
    for i in 1:N
        for j in 1+i:N
            if A[i, j] != 0
                τ[i, j] = τ_global
                τ[j, i] = τ_global
            end
        end
    end
    τ = round.(Integer, τ ./ Δt)

    # storing parameters
    save_ratio = 20
    no_saves = round(Integer, no_steps / save_ratio)
    store_θ = zeros(no_saves, N)
    θ_hist = rand(Uniform(0, 2 * π), maximum(τ) + 1, N)

    # init Kura object
    kura_sys = Kura_obj(σ, ω, A, θ_hist[end, :])

    # start simulation
    save_counter = 1
    for t in 1:no_steps

        coupling_sums = zeros(N)
                
        for i in 1:N

            coupling_sum = 0

            for j in 1:N
                coupling_sum += A[i, j] * sin(θ_hist[end - τ[i, j], j] - kura_sys.θ[i])
            end

            coupling_sums[i] = coupling_sum

        end

        # update θ short history
        θ_hist[1:end - 1, :] = θ_hist[2:end, :]
        θ_hist[end, :] = kura_sys.θ + (ω * Δt) + coupling_sums

        # update Kura obj
        kura_sys.θ = θ_hist[end, :]
        save_counter += 1

        # save θ
        if save_counter % save_ratio == 0
            store_θ[round(Integer, save_counter / save_ratio), :] = kura_sys.θ
        end
    end

    ####################################
    ############ save data #############

    params = HRGGnetwork[2]

    params["τ_global"] = τ_global
    params["nat_f"] = nat_f
    params["seedval"] = seedval
    params["Δt"] = Δt
    params["sim_time"] = sim_time

    # svae all results
    results = [store_θ, params, A]

    filepath = "/mnt/lustre/scratch/inf/ec627/data/HyperbolicNetworks/MetastabilitySoftHRGG/"

    fileseed = "Seed" * string(seedval)

    file_node = "_nodes_" * string(N)
    file_avg_degree = "_avg_degree_" * string(k)
    filecoupling = "_Coupling_K_" * string(K)
    file_τ_global = "_Delay_global_" * string(τ_global) * "_1e-3"
    file_γ = "_Gamma_" * string(γ)
    file_temp_ = "_Temp_" * string(T)

    filename = filepath * fileseed * file_node * file_avg_degree * filecoupling * file_τ_global * file_γ * file_temp_ * ".jld2"

    if test
        println("end test: ", filename)
    else
        save_object(filename, results)
    end
end
