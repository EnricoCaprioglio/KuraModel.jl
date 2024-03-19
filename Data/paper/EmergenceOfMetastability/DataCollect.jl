using Distributions
using Random
using JLD2
using LinearAlgebra

# The data I want to collect is:
# KOP of whole system: mean and standard deviation
# KOP of populations: mean and standard deviation
# KOP of modules at all timesteps

# 3D plot of how metastability changes as H and k change

# other parameters to save:
# H, k, K, n, B, β, sim_time, Δt, seedval

# define functions:
function SBMvar(n::AbstractArray, B::AbstractArray, H::Number, k::Number; K = 1)

	γ = (k + 1 - n[1]) / (n[1]*n[2] - n[1])

	return _SBMvar_constructor(n, B, H, γ; K = K)
	
end

function _partition_mat(n, B)

	if B[1] != B[2]*n[2]
		error("B₁ != B₂ × n₂, check model conditions")
	end
	if B[2] != B[3]*n[3]
		error("B₂ != B₃ × n₃, check model conditions")
	end

	N = prod(n)
	L = length(n)
	P = zeros(L - 1, N)
	
	for i in 1:(L-1)
		P[i, :] = vcat(
			[repeat([j], prod(n[1:i])) for j in 1:B[i]]...
		)
	end

	return P
	
end

function _SBMvar_constructor(n::AbstractArray, B::AbstractArray, H::Number, γ::Number; K = 1)

	if length(n) != length(B)
		error("vectors n and B mismatch")
	end
	if γ < 0 || γ > 1
		error("require 0 ≤ γ ≤ 1")
	end
	if B[1] != B[2]*n[2]
		error("B₁ != B₂ × n₂, check model conditions")
	end
	if B[2] != B[3]*n[3]
		error("B₂ != B₃ × n₃, check model conditions")
	end

	N = prod(n) # number of nodes
	L = length(n) # number of layers
	
	# store partition and adjaecncy matrix
	P = _partition_mat(n, B)
	A = zeros(N, N)

	# compute edge probabilities
	p₁ = (n[1]*H*γ) / (2*(n[1]-1)) + (1 - (n[1]*1*γ) / (2*(n[1]-1)))
	p₂ = (1/2 + H/2) * γ
	p₃ = (1/2 - H/2) * γ

	for i in 1:N
		for j in i+1:N

			# connect if in the same l = 1 block only
			if P[:, i] == P[:, j] && rand() < p₁
				A[i, j] = K
				A[j, i] = K
				
			# connect if in the same l = 2 block only
			elseif P[2, i] == P[2, j] && rand() < p₂
				A[i, j] = K
				A[j, i] = K

			# connect if in the same l = 3 block only
			elseif P[2, i] != P[2, j] && rand() < p₃
				A[i, j] = K
				A[j, i] = K
				
			end
		end
	end

	return A, P, [p₁, p₂, p₃]
end

function KS_sim(A::AbstractArray, P::AbstractArray, β::Real, K::Real, ω::AbstractArray, Δt::Real, sim_time::Real; seedval = 123, noiseQ = false, save_ratio = 1)

	N = size(A)[1]
	
	# simulations settings
	Random.seed!(seedval)
	
	# lag parameter
	α = π / 2 - β  # lag parameter in the equations
	
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
	if noiseQ
	    noise_scale = noise_scale
	else
	    noise_scale = 0
	end
	
	# simulation parameters
	steps = (0.0+Δt):Δt:sim_time
	no_steps = length(steps)
	
	# storing parameters
	no_saves = round(Integer, no_steps / save_ratio)
	θ_now = rand(Uniform(-π, π), N)  # random init conditions
	θs = zeros(no_saves, N)
	θs[1, :] = θ_now
	
	save_counter = 1
	
	for t in 2:no_steps
	    
	    # update phases
	    θj_θi_mat = (repeat(θ_now',N) - repeat(θ_now',N)') - α_mat
	    setindex!.(Ref(θj_θi_mat), 0.0, 1:N, 1:N) # set diagonal elements to zero 
	
	    k1 = map(sum, eachrow(A .* sin.(θj_θi_mat)))
	    θ_now += Δt .* (ω + k1) + noise_scale*(rand(Normal(0,1),N))*sqrt(Δt)
	    save_counter += 1
	
	    # save θ
	    if save_counter % save_ratio == 0
	        θs[round(Integer, save_counter / save_ratio), :] = θ_now
	    end
	    
	end
	
	params = Dict(
	    "α" => α,
	    "β" => β,
	    "K" => K,
	    # "A" => A,
		# "ω" => ω,
	    "seedval" => seedval,
	    "Δt" => Δt,
	    "sim_time" => sim_time,
	    "noiseQ" => noiseQ,
	)

	return θs, params
end

function macro_op(θs::AbstractMatrix)
	# num of cols = N
	N = length(θs[1, : ])
    map(θ -> abs(sum(exp.(im * θ) / N)), eachrow(θs))
end

# set folder path
folderpath = "/mnt/lustre/scratch/inf/ec627/data/HierarchicalChimera/paper_data/DataCollect/n1616/"

# job input
test = true
if test
    seedval = 123
else
    seedval = parse(Int, ARGS[1]) # this is the job input, setting the seedvalue
end
Random.seed!(seedval)

if test
    H_range = 0.4
else    
    H_range = 0.00:0.01:1
end
for H in H_range
    
    # hSBM parameters
    n = [16, 16, 2]
    B = [32, 2, 1]
    N = prod(n)
    k = N / 3
    K = 50

    # simulation parameters
    β = 0.1
    ω = repeat([1], N)
    save_ratio = 1
    if test
        sim_time = 1
        Δt = 1e-3
        relax = round(Integer, 0.1 / Δt)
    else
        sim_time = 55
        Δt = 1e-4
        relax = round(Integer, 5 / Δt)
    end
    steps = (0.0+Δt):Δt:sim_time
    no_steps = length(steps)

    # build network model
    Random.seed!(seedval); A, P, p = SBMvar(n, B, H, k; K = K)
    A = A ./ k # rescale coupling by average degree
    
    # simulation
    θs, params = KS_sim(A, P, β, K, ω, Δt, sim_time; seedval = seedval, noiseQ = false, save_ratio = 1)

    ## collect data
    save_ratio = 2
    # KOP whole system
    KOP_whole_mean = mean(macro_op(θs[relax:save_ratio:end, :]))
    KOP_whole_std = std(macro_op(θs[relax:save_ratio:end, :]))

    # KOP populations
    KOP_pop_1_mean = mean(macro_op(θs[relax:save_ratio:end, 1:n[1]*n[2]]))
    KOP_pop_1_std = std(macro_op(θs[relax:save_ratio:end, 1:n[1]*n[2]]))
    KOP_pop_2_mean = mean(macro_op(θs[relax:save_ratio:end, n[1]*n[2]+1:end]))
    KOP_pop_2_std = std(macro_op(θs[relax:save_ratio:end, n[1]*n[2]+1:end]))

    # KOP modules
    steps_to_save = length(θs[relax:save_ratio:end, 1])
    KOP_modules = zeros(steps_to_save, B[1])

    for i in 1:B[1]
        KOP_modules[:, i] = macro_op(θs[relax:save_ratio:end, n[1]*(i-1)+1:n[1]*i])
    end

    # other parameters to save:
    # H, k, K, n, B, β, sim_time, Δt, seedval
    params["k"] = k
    # params["n"] = n
    # params["B"] = B
    params["H"] = H
    params["save_ratio"] = save_ratio

    # save results
    results = [KOP_whole_mean, KOP_whole_std, KOP_pop_1_mean, KOP_pop_1_std, KOP_pop_2_mean, KOP_pop_2_std, KOP_modules, params]

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