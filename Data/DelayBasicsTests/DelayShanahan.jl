using Kuramodel
using Distributions
using Random
using JLD2

function get_partition(M, C)

	return vcat([repeat([i], M) for i in 1:C]...)
	
end

function get_pint(C, M, c_int)
    return c_int / (M - 1)
end

function get_pext(C, M, c_ext)
    return c_ext / (M * (C - 1))
end

function get_comm_mat(M, C, k, a, b, c_int)

    # init adj matrix and utilities
    local N = C * M
    A = zeros(N, N)
    local partition_vec = vcat([repeat([i], M) for i in 1:C]...)

    # couplings
    local μ = b * (a / k)
    local ν = (1 - b) * (a / k)
    
    # connection probabilities
    local c_ext = k - c_int
    local p_int = get_pint(C, M, c_int)
    local p_ext = get_pext(C, M, c_ext)
    
    for i in 1:N
        for j in i+1:N
            if partition_vec[i] == partition_vec[j] && rand() < p_int
                A[i, j] = μ
                A[j, i] = μ
            elseif partition_vec[i] != partition_vec[j] && rand() < p_ext
                A[i, j] = ν
                A[j, i] = ν
            end
        end
    end
    
    return A

end

seedval = 123
Random.seed!(seedval)

for τ_global in 3:1:7
    for c_int in [1,2,3,5,6,7]
        # network parameters
        M = 32
        C = 8
        N = M * C
        a = parse(Int, ARGS[1]) / 10
        b = 0.6
        k = 8
        partition_vec = get_partition(M, C)
        A = get_comm_mat(M, C, k, a, b, c_int)

        # oscillators parameters
        nat_f = 40  # Hz
        ω = repeat([2 * π * nat_f], N)
        σ = repeat([1], N)

        # simulation parameters
        Δt = 1e-4 # 0.1 ms
        A = A .* Δt  # rescaled by the time step already
        sim_time = 20  # seconds
        steps = (0.0 + Δt):Δt:sim_time
        no_steps = length(steps)

        # delays parameters
        # τ_global = 3
        τ_int = τ_global * 1e-3 # ms
        τ_ext = τ_global * 1e-3 # ms
        τ = zeros(N, N)
        for i in 1:N
            for j in i+1:N
                if A[i, j] != 0
                    if partition_vec[i] == partition_vec[j]
                        τ[i, j] = τ_int
                        τ[j, i] = τ_int
                    else
                        τ[i, j] = τ_ext
                        τ[j, i] = τ_ext
                    end
                end
            end
        end
        τ = round.(Integer, τ ./ Δt)

        # storing parameters
        save_ratio = 10
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

            # store_θ[t, :] = kura_sys.θ # this saves everything

        end


        ####################################
        ############ save data #############

        params = Dict(
            "M" => M,
            "C" => C,
            "N" => N,
            "a" => a,
            "b" => b,
            "k" => k,
            "c_int" => c_int,
            "τ_int" => τ_int,
            "τ_ext" => τ_ext,
            "nat_f" => nat_f,
            "seedval" => seedval,
            "Δt" => Δt,
            "sim_time" => sim_time,
            "A" => A
        )

        results = [store_θ, params]

        filepath = "/mnt/lustre/scratch/inf/ec627/data/DelayBasicsTests/output/"

        fileseed = "Seed" * string(seedval)

        filecoupling = "Coupling_a_" * string(a)
        file_extr_conn = "Ext_conn_" * string(c_int)
        file_τ_global = "Delay_global_" * string(τ_global) * "1e-3"

        filename = filepath * fileseed * filecoupling * file_extr_conn * file_τ_global * ".jld2"

        save_object(filename, results)
    end
end