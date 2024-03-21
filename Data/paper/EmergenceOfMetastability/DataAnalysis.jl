using Distributions
using Random
using JLD2
using LinearAlgebra

folderpath = "/mnt/lustre/scratch/inf/ec627/data/HierarchicalChimera/paper_data/DataCollect/n1608fast/"
cd(folderpath)
filenames = readdir()[1:end]

# Define a regular expression pattern to match the numbers
pattern = r"seed_(\d+)_beta_([\d.]+)_H_([\d.]+)_k_([\d.]+)\.jld2"

store_seeds = zeros(2, length(1:1:100))

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

        if k > 50
            store_seeds[1, seed] = seed
        else
            store_seeds[2, seed] = seed
        end

    end
end