using Distributions
using Random
using JLD2
using LinearAlgebra

folderpath = "/Users/ec627/Documents/Data/HierarchicalChimera/DataCollect/n1608fast/k_64"
H_range = collect(0.00:0.01:1)
no_seeds = 1

coal_entropy_store = zeros(no_seeds, lastindex(H_range))

cd(folderpath)
filenames = readdir()[1:end]

# Define a regular expression pattern to match the numbers
pattern = r"seed_(\d+)_beta_([\d.]+)_H_([\d.]+)_k_([\d.]+)\.jld2"

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

        if seed ≤ no_seeds

            KOP_whole_mean, KOP_whole_std, KOP_pop_1_mean, KOP_pop_1_std, KOP_pop_2_mean, KOP_pop_2_std, KOP_modules, params = load_object(filename)

            coalitions = get_coalitions(KOP_modules[1:2:end, :]'; γ = 0.9)
            coal_entropy_store[seed, round(Integer, H*100)+1] = coal_entropy(coalitions)
            
        end
    end
end

plt = plot(
    H_range,
    mean(coal_entropy_store, dims = 1)[1, :],
    xlabel = "H",
    ylabel = "Coalition Entropy",
    label = ""
)

savefig(plt, "/Users/ec627/Documents/Data/HierarchicalChimera/DataCollect/coalitions/CoalitionEntropies_γ_8.png")