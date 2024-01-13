using Kuramodel
using Distributions
using Random
using JLD2
using LinearAlgebra

path = "/mnt/lustre/scratch/inf/ec627/data/HierarchicalChimera/paper_data/raw_data"

test = true
if test
    seedval = 123
else
    seedval = parse(Int, ARGS[1])
end

Random.seed!(seedval)
	
# hSBM parameters
H = 0.4
n = [16, 8, 2]
B = [16, 2, 1]
N = prod(n)
k = N / 5
K = 1

# build network model
A, P, p = SBMvar(n, B, H, k; K = K)

# simulation parameters
β = 0.05
ω = repeat([1], N)
sim_time = 30
Δt = 1e-4
save_ratio = 1

# simulation
θs, params = KS_sim(A, P, β, K, ω, Δt, sim_time; seedval = seedval_KS_1, noiseQ = false, save_ratio = save_ratio)

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
relax = 5000
mean_pop1 = round(mean(macro_op(θs[relax+1:end, 1:Integer(N/2)])), digits = 3)
std_pop1 = round(std(macro_op(θs[relax+1:end, 1:Integer(N/2)])), digits = 3)
mean_pop2 = round(mean(macro_op(θs[relax+1:end, Integer(N/2):N])), digits = 3)
std_pop2 = round(std(macro_op(θs[relax+1:end, Integer(N/2):N])), digits = 3)

