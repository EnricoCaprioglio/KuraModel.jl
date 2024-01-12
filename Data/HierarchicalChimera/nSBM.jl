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

