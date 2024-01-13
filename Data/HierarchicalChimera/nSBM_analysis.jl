using Kuramodel
using JLD2

# set relevant paths
folderpath = "/mnt/lustre/scratch/inf/ec627/data/HierarchicalChimera/paper_data/raw_data"
filenames = readdir()[2:end]

# each file contains: results and params
# results = θs
# params = Dict(
# 	"n" => n,
# 	"B" => B,
# 	"N" => N,
# 	"α" => α,
# 	"β" => β,
# 	"K" => K,
# 	"B" => B,
# 	"A" => A,
# 	"seedval" => seedval,
# 	"Δt" => Δt,
# 	"sim_time" => sim_time,
# 	"H" => H,
# 	"P" => P
# )

# example filename
# seed_1_beta_0.05_H_0.45.jld2

# select params of interest
fileseed = [1,2,3]
filebeta = [0.05]
fileH = [0.45]

# the filename is then called using the following structure
# filename = folderpath * fileseed * filebeta * fileH * ".jld2"

for filename in filenames
    load_object(filename)
    local seed = parse(Int, filename[6])
    local β = parse(Int, filename[13:16])
    local H = parse(Int, filename[20:23])

    printl(filename)
    println("seed = $(seed), β = $β, H = $H")
end