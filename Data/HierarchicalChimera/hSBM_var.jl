using Kuramodel
using Distributions
using Random
using JLD2
using LinearAlgebra

# simulations settings
test = true
if test
    seedval = 123
else
    seedval = parse(Int, ARGS[1])
end

Random.seed!(seedval)

"""
	function SBMvar(n::AbstractArray, B::AbstractArray, H::Number; K = 1)

Function to crate the nested SBM variation with 3 hierarchical levels.

Inputs:

	n::AbstractArray 		vector of length 3, number of l-1 blocks in each l block
	B::AbstractArray 		number of blocks in at layer l
	H::Number 				value of hierarchical parameter

Outputs:

	A::AbstractMatrix 		 		adjacency matrix
	P::AbstractMatrix 				partition matrix
	[p₁, p₂, p₃]::AbstractArray 	probabilities of connection at layer l


"""
function SBMvar(n::AbstractArray, B::AbstractArray, H::Number; K = 1)

	if length(n) != length(B)
		error("vectors n and B mismatch")
	end

	N = prod(n) # number of nodes
	L = length(n) # number of layers
	
	# store partition and adjaecncy matrix
	P = zeros(L - 1, N)
	A = zeros(N, N)

	# compute edge probabilities
	ν = (1 - H) / 2
	μ = 1 - ν
	constant = 1 - (1 / (2*n[1] - 2))
	p₁ = (1 / (2*n[1] - 2)) * H + constant
	p₂ = μ / (K * n[1])
	p₃ = ν / (K * n[1])

	# create partition
	for i in 1:(L-1)
		P[i, :] = vcat(
			[repeat([j], prod(n[1:i])) for j in 1:B[i]]...
		)
	end

	for i in 1:N
		for j in i+1:N

			# connect if in the same module
			if P[:, i] == P[:, j] && rand() < p₁
				A[i, j] = K
				A[j, i] = K
				
			# connect with prob p_int if same population
			elseif P[2, i] == P[2, j] && rand() < p₂
				A[i, j] = K
				A[j, i] = K

			# connect with prob p_ext if diff population
			elseif P[2, i] != P[2, j] && rand() < p₃
				A[i, j] = K
				A[j, i] = K
				
			end
		end
	end

	return A, P, [p₁, p₂, p₃]
end

# hSBM parameters
H = 0.35
n = [64, 32, 2]
B = [64, 2, 1]
N = prod(n)
K = 1

# build network model
SBMvar_res = SBMvar(n, B, H; K = K)
A = SBMvar_res[1]
P = SBMvar_res[2]

# lag parameter
β = 0.1
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
ω = repeat([1], N) # identical
noiseQ = false# noise? (default = false)
if noiseQ
    noise_scale = 0.1
else
    noise_scale = 0
end

# simulation parameters
Δt = 1e-3
if test
    sim_time = .1
else
    sim_time = 30 # seconds
end
steps = (0.0+Δt):Δt:sim_time
no_steps = length(steps)

# storing parameters
save_ratio = 1
no_saves = round(Integer, no_steps / save_ratio)
global θ_now = rand(Uniform(-π, π), N)
θs = zeros(no_saves, N)
θs[1, :] = θ_now

global save_counter = 1

for t in 2:no_steps
    
    # update phases
    θj_θi_mat = (repeat(θ_now',N) - repeat(θ_now',N)') - α_mat
    setindex!.(Ref(θj_θi_mat), 0.0, 1:N, 1:N) # set diagonal elements to zero 

    k1 = map(sum, eachrow(A .* sin.(θj_θi_mat)))
    θ_now += Δt .* (ω + k1) # + noise_scale*(rand(Normal(0,1),N))*sqrt(Δt)
    save_counter += 1

    # save θ
    if save_counter % save_ratio == 0
        θs[round(Integer, save_counter / save_ratio), :] = θ_now
    end
    
end

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
    "noiseQ" => noiseQ,
    "P" => P
)

# save results
results = [θs, params]

folderpath = "/mnt/lustre/scratch/inf/ec627/data/HierarchicalChimera/hSBM_var/64_32_2_"

fileseed = "seed_" * string(seedval)
filebeta = "_beta_" * string(β)
fileH = "_H_" * string(H)

filename = folderpath * fileseed * filebeta * fileH

if test
    println("end test: ", filename * ".jld2")
    save_object(filename * "_TEST.jld2", results)
    println("File saved correctly: $(filename)")
else
    save_object(filename * ".jld2", results)
    println("File saved correctly: $(filename)")
end