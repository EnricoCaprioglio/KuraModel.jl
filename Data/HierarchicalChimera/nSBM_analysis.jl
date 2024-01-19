using Kuramodel
using Distributions
using Random
using JLD2
using LinearAlgebra

# each file contains: results and params
# results = θs
# params = Dict(
# 	"n" => n, "B" => B, "N" => N, "α" => α, "β" => β, "K" => K,
# 	"B" => B, "A" => A, "seedval" => seedval, "Δt" => Δt,
# 	"sim_time" => sim_time, "H" => H, "P" => P
# )

# example filename
# seed_9_beta_0.05_H_0.9_k_51.2.jld2
# 1234567890123456789012345678901234
# seed_3_beta_0.05_H_0.88_k_85.3.jld2
# 12345678901234567890123456789012345

# possible values
# H ∈ [0.00:0.01:1.00]
# k ∈ (51.2, 85.3)
# β ∈ (0.05)

# select params of interest
seed_range = parse(Int, ARGS[1])
β_range = 0.05; k_range = 85.3
H_range = collect(0.00:0.01:1.00)

# set relevant paths
folderpath = "/mnt/lustre/scratch/inf/ec627/data/HierarchicalChimera/paper_data/raw_data/k_" * string(k_range) * "/"
cd(folderpath)
filenames = readdir()[1:end]

# Define a regular expression pattern to match the numbers
pattern = r"seed_(\d+)_beta_([\d.]+)_H_([\d.]+)_k_([\d.]+)\.jld2"

# store data (each entry corresponds to a value of H, I added this since files are not ordered)
no_data = length(H_range)
means_pop_1 = zeros(no_data); means_pop_2 = zeros(no_data)
stds_pop_1 = zeros(no_data); stds_pop_2 = zeros(no_data)
means_modules_pop_1 = zeros(no_data); means_modules_pop_2 = zeros(no_data)
var_modules_pop_1 = zeros(no_data); var_modules_pop_2 = zeros(no_data)
Hs = zeros(no_data)

for filename in filenames
    # Match the pattern in the filename
    match_result = match(pattern, filename)

    # Extract values if there is a match
    if match_result !== nothing

        # Extract values from the matched groups
        seed = parse(Int, match_result[1])
        β = parse(Float64, match_result[2])
        H = parse(Float64, match_result[3])
        k = parse(Float64, match_result[4])

        if seed in seed_range && β == β_range && k == k_range

            # add to list of H studied
            append!(Hs, H)
            id = findfirst(x -> x == H, H_range)

            # load data
            θs, params = load_object(folderpath * filename)
            B = params["B"]; n = params["n"]; N = prod(n)
            relax = round(Integer, (1 / params["Δt"]) * 5)

            # layer 2 analysis (populations)
            pop_1_end = round(Integer, N/2)
            means_pop_1[id] = mean(macro_op(θs[relax+1:end, 1:pop_1_end]))
            means_pop_2[id] = mean(macro_op(θs[relax+1:end, pop_1_end+1:end]))
            stds_pop_1[id] = std(macro_op(θs[relax+1:end, 1:pop_1_end]))
            stds_pop_2[id] = std(macro_op(θs[relax+1:end, pop_1_end+1:end]))

            # layer 1 analysis (modules)
            modules_KOP_means = zeros(B[1]) # store mean KOP each module
            modules_KOP_std = zeros(B[1]) # store var KOP each module

            for i in 1:B[1] # this loops through each module (regardless of population)
                modules_macro_op = macro_op(θs[relax+1:end, n[1]*(i-1)+1:n[1]*i])
                modules_KOP_means[i] = mean(modules_macro_op)
                modules_KOP_std[i] = std(modules_macro_op)
            end

            modules_pop_1_end = Integer(B[1]/2)
            # get average of all modules in a single population
            means_modules_pop_1[id] = mean(modules_KOP_means[1:modules_pop_1_end])
            means_modules_pop_2[id] = mean(modules_KOP_means[modules_pop_1_end+1:end])
            var_modules_pop_1[id] = mean(modules_KOP_std[1:modules_pop_1_end])^2
            var_modules_pop_2[id] = mean(modules_KOP_std[modules_pop_1_end+1:end])^2
        end

        # Print or store the extracted values
        # println("seed = $seed, beta = $beta, H = $H, k = $k")
    else
        println("No match found in the filename: $(filename)")
    end
end

# save data
savefolder_path = "/mnt/lustre/scratch/inf/ec627/data/HierarchicalChimera/paper_data/data_analysis/"

results = Dict(
    "means_pop_1" => means_pop_1,
    "means_pop_2" => means_pop_2,
    "stds_pop_1" => stds_pop_1,
    "stds_pop_2" => stds_pop_2,
    "means_modules_pop_1" => means_modules_pop_1,
    "means_modules_pop_2" => means_modules_pop_2,
    "var_modules_pop_1" => var_modules_pop_1,
    "var_modules_pop_2" => var_modules_pop_2,
    "H" => Hs
)

savefilename = "_beta_" * string(β_range) * "_k_" * string(k_range) * ".jld2"

save_object(savefolder_path * savefilename, results)

println("Results saved at: $(savefolder_path * savefilename)")