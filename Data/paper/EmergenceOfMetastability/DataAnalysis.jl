using Distributions
using Random
using JLD2
using LinearAlgebra

########################################################
# change the path and filename if k_desired is changed #
########################################################

folderpath = "/mnt/lustre/scratch/inf/ec627/data/HierarchicalChimera/paper_data/DataCollect/n1608fast/k_42/"
k_desired = 42 # 51, 42, 64
no_seeds = 100
H_range = collect(0.00:0.01:1)
n₂ = 8
test = false

	
cd(folderpath)
if test
    filenames = readdir()[1:10]
else
    filenames = readdir()[1:end]
end

# Define a regular expression pattern to match the numbers
# pattern = r"seed_(\d+)_beta_([\d.]+)_H_([\d.]+)_k_([\d.]+)\.jld2"

whole_mean_store = zeros(no_seeds, lastindex(H_range))
whole_std_store = zeros(no_seeds, lastindex(H_range))
pop_1_mean_store = zeros(no_seeds, lastindex(H_range))
pop_1_std_store = zeros(no_seeds, lastindex(H_range))
pop_2_mean_store = zeros(no_seeds, lastindex(H_range))
pop_2_std_store = zeros(no_seeds, lastindex(H_range))
modules_means_pop_1_store = zeros(no_seeds, lastindex(H_range))
modules_means_pop_2_store = zeros(no_seeds, lastindex(H_range))
metastability_pop_1_store = zeros(no_seeds, lastindex(H_range))
metastability_pop_2_store = zeros(no_seeds, lastindex(H_range))

println(length(filenames))

for filename in filenames

# # Match the pattern in the filename
# match_result = match(pattern, filename)

# # Extract values if there is a match
# if match_result !== nothing

#     # Extract values from the matched groups
#     seed = parse(Int, match_result[1])
#     β = parse(Float64, match_result[2])
#     H = parse(Float64, match_result[3])
#     k = parse(Float64, match_result[4])

#     if k == k_desired
        
    KOP_whole_mean, KOP_whole_std, KOP_pop_1_mean, KOP_pop_1_std, KOP_pop_2_mean, KOP_pop_2_std, KOP_modules, params = load_object(folderpath * filename)

    H = params["H"]
    seed = params["seedval"]

    # compute the mean KOP between modules in pop 1 or 2
    pop1_modules_mean_KOP = mean(mean.([KOP_modules[:, i] for i in 1:n₂]))
    pop1_metastability = mean(std.([KOP_modules[:, i] for i in 1:n₂]))
    pop2_modules_mean_KOP = mean(mean.([KOP_modules[:, i] for i in n₂+1:n₂*2]))
    pop2_metastability = mean(std.([KOP_modules[:, i] for i in n₂+1:n₂*2]))

    ## store data
    # whole system store
    whole_mean_store[seed, round(Integer, H*100)+1] = KOP_whole_mean
    whole_std_store[seed, round(Integer, H*100)+1] = KOP_whole_std
    
    # populations store
    pop_1_mean_store[seed, round(Integer, H*100)+1] = KOP_pop_1_mean
    pop_1_std_store[seed, round(Integer, H*100)+1] = KOP_pop_1_std
    pop_2_mean_store[seed, round(Integer, H*100)+1] = KOP_pop_2_mean
    pop_2_std_store[seed, round(Integer, H*100)+1] = KOP_pop_2_std

    # modules store
    modules_means_pop_1_store[seed, round(Integer, H*100)+1] = pop1_modules_mean_KOP
    metastability_pop_1_store[seed, round(Integer, H*100)+1] = pop1_metastability
    modules_means_pop_2_store[seed, round(Integer, H*100)+1] = pop2_modules_mean_KOP
    metastability_pop_2_store[seed, round(Integer, H*100)+1] = pop2_metastability

end

## stable and unstable populations
# stable and unstable populations
stable_mean = zeros(size(pop_1_mean_store))
stable_std = zeros(size(pop_1_std_store))
unstable_mean = zeros(size(pop_1_mean_store))
unstable_std = zeros(size(pop_1_std_store))

# stable and unstable modules
stable_modules_mean = zeros(size(modules_means_pop_1_store))
stable_modules_metastability = zeros(size(modules_means_pop_1_store))
unstable_modules_mean = zeros(size(metastability_pop_1_store))
unstable_modules_metastability = zeros(size(metastability_pop_1_store))

# assign stable and unstable labels to populations
for j in 1:no_seeds
    for i in 1:lastindex(H_range)
        
        δ = 0.1

        # assign stable and unstable populations
        if pop_1_mean_store[j, i] > pop_2_mean_store[j, i] # pop 1 stable
            
            stable_mean[j, i] = pop_1_mean_store[j, i]
            stable_std[j, i] = pop_1_std_store[j, i] 
            stable_modules_mean[j, i] = modules_means_pop_1_store[j, i] 
            stable_modules_metastability[j, i] = metastability_pop_1_store[j, i] 

            unstable_mean[j, i] = pop_2_mean_store[j, i]
            unstable_std[j, i] = pop_2_std_store[j, i]
            unstable_modules_mean[j, i] = modules_means_pop_2_store[j, i]
            unstable_modules_metastability[j, i] = metastability_pop_2_store[j, i]
            
        else

            stable_mean[j, i] = pop_2_mean_store[j, i]
            stable_std[j, i] = pop_2_std_store[j, i] 
            stable_modules_mean[j, i] = modules_means_pop_2_store[j, i] 
            stable_modules_metastability[j, i] = metastability_pop_2_store[j, i] 

            unstable_mean[j, i] = pop_1_mean_store[j, i]
            unstable_std[j, i] = pop_1_std_store[j, i]
            unstable_modules_mean[j, i] = modules_means_pop_1_store[j, i]
            unstable_modules_metastability[j, i] = metastability_pop_1_store[j, i]
            
        end
    end
end

# population metastability index
pop_metastability = zeros(length(stable_std[:, 1])*2, length(stable_std[1, :]))
pop_metastability[1:length(stable_std[:, 1]), :] = stable_std
pop_metastability[length(stable_std[:, 1])+1:end, :] = unstable_std

# modules metatsability index
modules_metastability = zeros(length(stable_modules_metastability[:, 1])*2, length(stable_modules_metastability[1, :]))
modules_metastability[1:length(stable_std[:, 1]), :] = stable_modules_metastability
modules_metastability[length(stable_std[:, 1])+1:end, :] = unstable_modules_metastability


data_to_plot = Dict(
    "mean_whole_mean_store" => mean(whole_mean_store, dims = 1)[1, :],
    "std_whole_mean_store" => std(whole_mean_store, dims = 1)[1, :],
    "mean_stable_mean" => mean(stable_mean, dims = 1)[1, :],
    "std_stable_mean" => mean(stable_mean, dims = 1)[1, :],
    "mean_unstable_mean" => mean(unstable_mean, dims = 1)[1, :],
    "std_unstable_mean" => std(unstable_mean, dims = 1)[1, :],
    "mean_stable_std" => mean(stable_std, dims = 1)[1, :],
    "std_stable_std" => std(stable_std, dims = 1)[1, :],
    "mean_unstable_std" => mean(unstable_std, dims = 1)[1, :],
    "std_unstable_std" => std(unstable_std, dims = 1)[1, :],
    "mean_pop_metastability" => mean(pop_metastability, dims = 1)[1, :],
    "std_pop_metastability" => std(pop_metastability, dims = 1)[1, :],
    "mean_stable_modules_mean" => mean(stable_modules_mean, dims = 1)[1, :],
    "std_stable_modules_mean" => std(stable_modules_mean, dims = 1)[1, :],
    "mean_unstable_modules_mean" => mean(unstable_modules_mean, dims = 1)[1, :],
    "std_unstable_modules_mean" => std(unstable_modules_mean, dims = 1)[1, :],
    "mean_stable_modules_metastability" => mean(stable_modules_metastability .^2, dims = 1)[1, :],
    "std_stable_modules_metastability" => std(stable_modules_metastability .^2, dims = 1)[1, :],
    "mean_unstable_modules_metastability" => mean(unstable_modules_metastability .^2, dims = 1)[1, :],
    "std_unstable_modules_metastability" => std(unstable_modules_metastability .^2, dims = 1)[1, :],
    "whole_modules_mean_metastability" => mean(modules_metastability .^2, dims = 1)[1, :],
    "whole_modules_std_metastability" => std(modules_metastability .^2, dims = 1)[1, :],
)


filename = "_no_seeds_" * string(no_seeds) * "_k_" * string(k_desired) * ".jld2"
if test
    save_object("/mnt/lustre/scratch/inf/ec627/data/HierarchicalChimera/paper_data/DataCollect/n1608fast/k_42_plots/TEST" * filename, data_to_plot)
    println("File saved correctly: $(filename)")
else
    save_object("/mnt/lustre/scratch/inf/ec627/data/HierarchicalChimera/paper_data/DataCollect/n1608fast/k_42_plots/" * filename, data_to_plot)
    println("File saved correctly: $(filename)")
end