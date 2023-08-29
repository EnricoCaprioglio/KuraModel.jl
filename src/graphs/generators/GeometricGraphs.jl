using Kuramodel

"""
	function hyperbolic_graph(graphtype, N, k, γ, T; ζ = 1, ϵ = 0.01, max_iter = 100)

Function used to generate Hyperbolic Random Geometric Graphs.

The input parameters are:
- the type of graph: graphtype = :HRGG or :softHRGG
- number of nodes `N`;
- target average degree `k`;
- target degree distribution esponent `γ`;
- temperature `T`.

Default parameters:
- curvature `K = -ζ ² = - 1`;
- precision `ϵ` to calculate the radius `R` of the Pioncare disc;
- maximum number of iterations `max_iter` to calculate `R`.

Output  a vector of length 4 containing:
 - `[1]` = the adjacency matrix `A` of the generated network ;
- `[2]` = a dictionary collecting the parameters used to generate the network;
- `[3]` = radial coordinates `r_i`
- `[4]` = angular coordinates `θ_i`
"""
function hyperbolic_graph(graphtype, N, k, γ, T; ζ = 1, ϵ = 0.01, max_iter = 100)

	if γ < 2
		error("γ > 3 is required")
	end
	
	if graphtype == :HRGG
		if T != 0.0
			error("To generate HRGG you need T -> 0")
		end
	end

	if graphtype == :softHRGG
		if T == 0.0
			error("To generate a soft HRGG you need 0 < T < ∞")
		end
	end

	# store parameters here
	params = Dict(
		"k" => k,
		"γ" => γ,
		"T" => T,
		"N" => N,
		"ζ" => ζ,
		"ϕ1" => 0.,
		"ϵ" => ϵ,
		"max_iter" => max_iter
	)

	# step 1: find α
	α = get_alpha(T, γ; ζ = ζ)
	params["α"] = α
	
	# step 2: find R
	ϕ1 = params["ϕ1"]
	if graphtype == :HRGG
		R = get_R_HRGG(k, ϵ, N, α, ζ, T, ϕ1; R1 = 0.001, R2 = 15, max_iter = max_iter)
		params["R"] = R
	elseif graphtype == :softHRGG
		R = get_R_softHRGG(k, ϵ, N, α, ζ, T, ϕ1; R1 = 0.001, R2 = 15, max_iter = max_iter)
		params["R"] = R
	end

	# step 3: sample coords
	# radial coords
	r_i = sample_r(R, N, α)
    # angular coords
	ϕ_i = sort(sample_ϕ(N; distr = :uniform))

	# step 4: connect
	if graphtype == :HRGG
		A = connect_nodes_HRGG(r_i, ϕ_i, R; ζ = ζ)
	elseif graphtype == :softHRGG
		A = connect_nodes_softHRGG(r_i, ϕ_i, R, T; ζ = 1)
	end

	return [A, params, r_i, ϕ_i]
end