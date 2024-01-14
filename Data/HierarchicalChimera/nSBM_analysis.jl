using Kuramodel
using Distributions
using Random
using JLD2
using LinearAlgebra

# set relevant paths
folderpath = "/mnt/lustre/scratch/inf/ec627/data/HierarchicalChimera/paper_data/raw_data"
cd(folderpath)
filenames = readdir()[1:end]

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

# select params of interest
seed_range = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
β_range = 0.05; H_range = 0.45; k_range = 51.2

# Define a regular expression pattern to match the numbers
pattern = r"seed_(\d+)_beta_([\d.]+)_H_([\d.]+)_k_([\d.]+)\.jld2"

# store data
means_pop_1 = []; means_pop_2 = []
stds_pop_1 = []; stds_pop_2 = []
means_modules_pop_1 = []; means_modules_pop_2 = []
var_modules_pop_1 = []; var_modules_pop_2 = []

for filename in filenames
    # Match the pattern in the filename
    match_result = match(pattern, filename)

    # Extract values if there is a match
    if match_result !== nothing

        # Extract values from the matched groups
        seed = parse(Int, match_result[1])
        beta = parse(Float64, match_result[2])
        H = parse(Float64, match_result[3])
        k = parse(Float64, match_result[4])

        if seed in seed_range && β == β_range && H == H_range && k == k_range

            # load data
            θs, params = load_object(folderpath * filename)
            B = params["B"]; n = params["n"]; N = prod(N)
            relax = round(Integer, (1 / params["Δt"]) * 5)

            # layer 2 analysis (populations)
            append!(means_pop_1, mean(macro_op(θs[relax+1:end, 1:round(Integer, N/2)])))
            append!(means_pop_2, mean(macro_op(θs[relax+1:end, round(Integer, N/2)+1:end])))
            append!(stds_pop_1, std(macro_op(θs[relax+1:end, 1:round(Integer, N/2)])))
            append!(stds_pop_2, std(macro_op(θs[relax+1:end, round(Integer, N/2)+1:end])))

            # layer 1 analysis (modules)
            modules_KOP = zeros(B[1])
            modules_KOP_std = zeros(B[1])
            modules_macro_op = zeros(length(θs[relax:end, 1]), B[1])

            for i in 1:B[1]
                modules_macro_op[:, i] = macro_op(θs[relax:end, n[1]*(i-1)+1:n[1]*i])
                modules_KOP[i] = mean(modules_macro_op[:, i])
                modules_KOP_std[i] = std(modules_macro_op[:, i])
            end

            append!(means_modules_pop_1, mean(modules_KOP[1:Integer(B[1]/2)]))
            append!(means_modules_pop_2, mean(modules_KOP[Integer(B[1]/2)+1:end]))
            append!(var_modules_pop_1, mean(modules_KOP_std[1:Integer(B[1]/2)])^2)
            append!(var_modules_pop_2, mean(modules_KOP_std[Integer(B[1]/2)+1:end])^2)
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
    "var_modules_pop_2" => var_modules_pop_2
)

savefilename = "_beta_" * string(β_range) * "_H_" * string(H_range) * "_k_" * string(k_range) * ".jld2"

save_object(savefolder_path * savefilename)