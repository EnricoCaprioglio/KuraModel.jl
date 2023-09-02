using Kuramodel
using JLD2
using Random
using Distributions
using LinearAlgebra
using MCIntegration

seedval = 123
Random.seed!(seedval)

function test_numerical_k(graphtype, N, γ, T, var_k, iter_val; ζ = 1, ϵ = 0.01, max_iter = 100)

	store_numerical_k = zeros(iter_val, length(var_k))
	
	for (iter, k) in enumerate(var_k)
	
		for i in 1:iter_val

			if graphtype == :HRGG
				
				HRGG_test = hyperbolic_graph(:HRGG, N, k, γ, 0.0; ζ = ζ, ϵ = ϵ, max_iter = max_iter)
				
			elseif graphtype == :softHRGG
				
				HRGG_test = hyperbolic_graph(:softHRGG, N, k, γ, T; ζ = ζ, ϵ = ϵ, max_iter = max_iter)
				
			end
			
			numerical_k = mean([sum(HRGG_test[1][i, :]) for i in 1:HRGG_test[2]["N"]])
			store_numerical_k[i, iter] = numerical_k
		end
	end
	
    out1 = var_k
    out2 = [mean(store_numerical_k[:, i]) for i in 1:length(var_k)]
    out3 = [std(store_numerical_k[:, i]) for i in 1:length(var_k)]

    params = Dict(
        "graphtype" => graphtype,
		"γ" => γ,
		"T" => T,
		"N" => N,
		"ζ" => ζ,
		"ϵ" => ϵ,
		"max_iter" => max_iter,
        "iter_val" => iter_val
	)
    
    return [out1, out2, out3, params]

	# println(
	# 	"Target average degrees: ", var_k, "\n",
	# 	"Numerical degrees obtained: ", [mean(store_numerical_k[:, i]) for i in 1:length(var_k)],
	# 	"\nNumerical degrees stds: ", [std(store_numerical_k[:, i]) for i in 1:length(var_k)]
	# )
	
end

test = false

######### parameters #########

# **Convert array job input to a number**
# taskID should be between 1 and 10
if test
	N = 100
else
	N = 100 * parse(Int, ARGS[1])  # number of nodes
end

graphtype = :softHRGG

if test
    var_k = [10, 20]
else
    var_k = collect(5:5:100)
end

γ = 2.5

if graphtype == :HRGG
    T = 0.0
elseif graphtype == :softHRGG
    T = 0.5
end

if test
    iter_val = 2
else
    iter_val = 100
end

ϵ = 0.1

######### computation #########

results = test_numerical_k(graphtype, N, γ, T, var_k, iter_val; ζ = 1, ϵ = ϵ, max_iter = 100)

####################################
############ save data #############

filepath = "/its/home/ec627/data/HyperbolicNetworks/DegreeTests/outputs/"

fileseed = "Seed" * string(seedval)
fileN = "noNodes" * string(N)
fileVarK = "varK5to100"
filegraphtype = "GraphHRGG"

# all other parameters are stored in results[3] i.e., the dictionary
# Dict(
#     "graphtype" => graphtype,
#     "γ" => γ,
#     "T" => T,
#     "N" => N,
#     "ζ" => ζ,
#     "ϵ" => ϵ,
#     "max_iter" => max_iter,
#     "iter_val" => iter_val
# )

if test
    filename = filepath * fileseed * fileN * fileVarK * filegraphtype * "TEST.jld2"
else
    filename = filepath * fileseed * fileN * fileVarK * filegraphtype * ".jld2"
end

println("Job finished, saving object.\nParameters used were:\n", results)

save_object(filename, results)

####################################
