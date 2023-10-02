using Kuramodel
using Distributions
using Random
using JLD2
using LinearAlgebra

test = true
noiseQ = false

if test
    job_ID = 3
else
    job_ID = parse(Int, ARGS[1])
end

# set seed value
seedval = job_ID
Random.seed!(seedval)

B = 0.25
ν = (1 - B) / 2
μ = 1 - ν
M = 32
C = 64
N = M * C * 2
K = 1

p_int = μ / (K * M)
p_ext = ν / (K * M)

c_int = p_int * M * (C - 1)
c_ext = p_ext * M * C

# adjacency matrix construction
partition = zeros(2, N)
partition[1, 1:(M*C)] = vcat([repeat([i], M) for i in 1:C]...)
partition[1, (M*C)+1:N] = vcat([repeat([i], M) for i in 1:C]...)
partition[2, :] = vcat([repeat([i], (M*C)) for i in 1:2]...)

# construct adjacency matrix
A = zeros(N, N)
for i in 1:N
    for j in i+1:N

        # connect if in the same module
        if partition[:, i] == partition[:, j]
            A[i, j] = K * 2
            A[j, i] = K * 2

        # connect with prob p_int if same population
        elseif partition[2, i] == partition[2, j] && rand() < p_int
            A[i, j] = K
            A[j, i] = K

        # connect with prob p_ext if diff population
        elseif partition[2, i] != partition[2, j] && rand() < p_ext
            A[i, j] = K
            A[j, i] = K
            
        end
    end
end

# lag parameters
β = 0.1
α = π / 2 - β  # lag parameter in the equations
α_mat = zeros(N, N)

# construct lag matrix
for i in 1:N
    for j in i+1:N

        # add lag if different module
        if A[i,j] != 0 && partition[:, i] != partition[:, j]
            α_mat[i, j] = α
            α_mat[j, i] = α
            
        end
    end
end

# oscillators parameters
ω = repeat([1], N)

# simulation parameters
Δt = 1e-3
sim_time = 20
steps = (0.0+Δt):Δt:sim_time
if test
    no_steps = 10
else
    no_steps = length(steps)
end
if noiseQ
    noise_scale = 0.1
else
    noise_scale = 0
end

# storing parameters
if test
    save_ratio = 1
else
    save_ratio = 10
end
no_saves = round(Integer, no_steps / save_ratio)
store_θ = zeros(no_saves, N)
global θ_now = rand(Uniform(-π, π), N)

# start simulation
global save_counter = 1
for t in 1:no_steps
    
    # update phases
    θj_θi_mat = (repeat(θ_now',N) - repeat(θ_now',N)') - α_mat
    setindex!.(Ref(θj_θi_mat), 0.0, 1:N, 1:N) # set diagonal elements to zero 

    k1 = map(sum, eachrow(A .* sin.(θj_θi_mat)))
    global θ_now += Δt .* (ω + k1) # + noise_scale*(rand(Normal(0,1),N))*sqrt(Δt)

    # save θ
    if save_counter % save_ratio == 0
        store_θ[round(Integer, save_counter / save_ratio), :] = θ_now
    end
    
    global save_counter += 1

end

params = Dict(
    "M" => M,
    "C" => C,
    "N" => N,
    "α" => α,
    "β" => β,
    "K" => K,
    "B" => B,
    "A" => A,
    "seedval" => seedval,
    "Δt" => Δt,
    "sim_time" => sim_time,
    "μ" => μ,
    "ν" => ν,
    "noiseQ" => noiseQ,

)

results = [store_θ, params]

folderpath = "/mnt/lustre/scratch/inf/ec627/data/HierarchicalChimera/VarM/"

fileseed = "Seed" * string(seedval)
filebeta = "_beta_" * string(β)
fileB = "_B_" * string(B)
fileM = "_M_" * string(M)
fileC = "_C_" * string(C)

filename = folderpath * fileseed * filebeta * fileB * fileM * fileC

if test
    println("end test: ", filename)
else
    save_object(filename, results)
end