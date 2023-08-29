using MCIntegration
using Distributions
using Random

"""
	function get_alpha(T, γ; ζ = 1)

Function to calculate α given:
- temperature T
- target exponent γ.

Curvature is assumed to be K = -ζ ² = -1.
Similarly, if ζ, T -> ∞ then η = 1 by default (configuration model).

"""
function get_alpha(T::Number, γ::Number; ζ = 1, η = 1)
	
	if γ < 2
		error("Gamma is smaller than 2")
	end
	
	if T == 1
		error("Temperture is exactly 1, critical point")
	end
		
	if T < 1
		return ζ * (γ - 1) / 2
	elseif T > 1
		return ζ * (γ - 1) / (2 * T)
    elseif T > 10000
        return η * (γ - 1) / 2
	end

end

"""
	function get_ξ(α::Number, ζ::Number)

Simple function to calculate ξ = (α/ζ)/(α/ζ - 1/2).
"""
function get_ξ(α::Number, ζ::Number)

	return (α/ζ)/(α/ζ - 1/2)

end

"""
    function px_HRGG(R::Number, x::Number)

Connection probability for Hyperbolic RGG with T -> 0.
"""
function px_HRGG(R::Number, x::Number)
	if R - x > 0
		return 1
	else
		return 0
	end
end

"""
    function px_softHRGG(R::Number, x::Number, T::Number, ζ::Number)

Connection probability for soft Hyperbolic RGG with T ∈ (0, ∞).
"""
function px_softHRGG(R::Number, x::Number, T::Number, ζ::Number)

	β = 1/T
	k1 = exp(β * (ζ / 2) * (x - R))
	
	return 1/ (1 + k1)
end

"""
    function integrand_HRGG(vars, config)

Integrand function used to solve numerically the integral to find the average degree of Hyperbolic RGG.
"""
function integrand_HRGG(vars, config)

	# variables to integrate
	r1, r2, ϕ2 = vars

	# parameters
	parameters = config.userdata
	
	T = parameters["T"]
	β = 1/T
	ζ = parameters["ζ"]
	α = parameters["α"]
	N = parameters["N"]
	ϕ1 = parameters["ϕ1"]
	R = parameters["R"]

	# eval hyper distance
	Δϕ = π - abs(π - abs(ϕ1 - ϕ2[1]))
	k1 = cosh(ζ * r1[1]) * cosh(ζ * r2[1])
	k2 = sinh(ζ * r1[1]) * sinh(ζ * r2[1]) * cos(Δϕ)
	x = acosh(k1 - k2) / ζ

	# eval ρ(r1) * ρ(r2) * p(x) = k3 * k2 * k1
	k₁ = px_HRGG(R, x)
	k₂ = α * (sinh(α * r2[1])) / (cosh(α * R) - 1)
	k₃ = α * (sinh(α * r1[1])) / (cosh(α * R) - 1)

	return (N/π) * k₃ * k₂ * k₁
end

"""
    function integrand_softHRGG(vars, config)

Integrand function used to solve numerically the integral to find the average degree of a soft Hyperbolic RGG.
"""
function integrand_softHRGG(vars, config)

	# variables to integrate
	r1, r2, ϕ2 = vars

	# parameters
	parameters = config.userdata
	
	T = parameters["T"]
	β = 1/T
	ζ = parameters["ζ"]
	α = parameters["α"]
	N = parameters["N"]
	ϕ1 = parameters["ϕ1"]
	R = parameters["R"]

	# eval hyper distance
	Δϕ = π - abs(π - abs(ϕ1 - ϕ2[1]))
	k1 = cosh(ζ * r1[1]) * cosh(ζ * r2[1])
	k2 = sinh(ζ * r1[1]) * sinh(ζ * r2[1]) * cos(Δϕ)
	x = acosh(k1 - k2) / ζ

	# eval ρ(r1) * ρ(r2) * p(x) = k3 * k2 * k1
	k₁ = px_softHRGG(R, x, T, ζ)
	k₂ = α * (sinh(α * r2[1])) / (cosh(α * R) - 1)
	k₃ = α * (sinh(α * r1[1])) / (cosh(α * R) - 1)

	return (N/π) * k₃ * k₂ * k₁
end

"""
	function get_avg_k_HRGG(R::Number, params)

Function to calculate the average degree k of the network.
It uses the package `MCIntegration.jl`

Input parameters:

- `R`, the radius of the Poincare' disc;
- `params`, a dictionary of parameters.

Note, `params` needs to contain the keys:
- `[ "T", "ζ", "α", "N", "ϕ1" ]`

Outputs:
- a vector `[avg_k, avg_k_std]` containing the integration result `avg_k` and the standard deviation `avg_k_std`.
"""
function get_avg_k_HRGG(R::Number, params)

	# change R if it's already in the parameters
	params["R"] = R
	
	# set variables
	r1r2ϕ2 = MCIntegration.Continuous(
		[
			(0.0, params["R"]),
			(0.0, params["R"]),
			(0.0, π)
		]
	)

	# integrate
	res_avg_k = integrate(
		integrand_HRGG;
		var = r1r2ϕ2,
		userdata = params
	)

	avg_k = res_avg_k.mean[1]
	avg_k_std = res_avg_k.stdev[1]

	return [avg_k, avg_k_std]
end

"""
	function get_avg_k_softHRGG(R::Number, params)

Function to calculate the average degree k of the network.
It uses the package `MCIntegration.jl`

Input parameters:

- `R`, the radius of the Poincare' disc;
- `params`, a dictionary of parameters.

Note, `params` needs to contain the keys:
- `[ "T", "ζ", "α", "N", "ϕ1" ]`

Outputs:
- a vector `[avg_k, avg_k_std]` containing the integration result `avg_k` and the standard deviation `avg_k_std`.
"""
function get_avg_k_softHRGG(R::Number, params)

	# change R if it's already in the parameters
	params["R"] = R
	
	# set variables
	r1r2ϕ2 = MCIntegration.Continuous(
		[
			(0.0, params["R"]),
			(0.0, params["R"]),
			(0.0, π)
		]
	)

	# integrate
	res_avg_k = integrate(
		integrand_softHRGG;
		var = r1r2ϕ2,
		userdata = params
	)

	avg_k = res_avg_k.mean[1]
	avg_k_std = res_avg_k.stdev[1]

	return [avg_k, avg_k_std]
end

"""
	function get_R_HRGG(k, ϵ, N, α, ζ, T, ϕ1; R1 = 0.001, R2 = 15, max_iter = 100)

Bisection method to find the required `R` given the target average degree `k`.

"""
function get_R_HRGG(k, ϵ, N, α, ζ, T, ϕ1; R1 = 0.001, R2 = 15, max_iter = 100)

	# store parameters
	params = Dict(
		"ζ" => ζ,
		"α" => α,
		"N" => N,
		"ϕ1" => ϕ1, 
		"T" => T
	)

	iter = 0
	k_temp = 10000
	# ϵ_collect = []
	
	while abs(k_temp - k) > ϵ && iter < max_iter

		# set mid point R
		params["R"] = (R1 + R2)/2
		R = params["R"]
		
		k_temp = get_avg_k_HRGG(R, params)[1]

		# update R1 and R2 accordingly
		if k_temp < k
			R1 = R1
			R2 = R
		else
			R1 = R
			R2 = R2
		end
		
		iter += 1
		# append!(ϵ_collect, abs(k_temp - k))
	end

	return (R1 + R2)/2 # iter, ϵ_collect
end

"""
	function get_R_softHRGG(k, ϵ, N, α, ζ, T, ϕ1; R1 = 0.001, R2 = 15, max_iter = 100)

Bisection method to find the required `R` given the target average degree `k`.

"""
function get_R_softHRGG(k, ϵ, N, α, ζ, T, ϕ1; R1 = 0.001, R2 = 15, max_iter = 100)

	# store parameters
	params = Dict(
		"ζ" => ζ,
		"α" => α,
		"N" => N,
		"ϕ1" => ϕ1, 
		"T" => T
	)

	iter = 0
	k_temp = 10000
	# ϵ_collect = []
	
	while abs(k_temp - k) > ϵ && iter < max_iter

		# set mid point R
		params["R"] = (R1 + R2)/2
		R = params["R"]
		
		k_temp = get_avg_k_softHRGG(R, params)[1]

		# update R1 and R2 accordingly
		if k_temp < k
			R1 = R1
			R2 = R
		else
			R1 = R
			R2 = R2
		end
		
		iter += 1
		# append!(ϵ_collect, abs(k_temp - k))
	end

	return (R1 + R2)/2 # iter, ϵ_collect
end

"""
	function sample_r(R::Number, N::Integer, α::Number)

Function used to sample radial coordinates rᵢ.
"""
function sample_r(R::Number, N::Integer, α::Number)
	
	Ui = rand(N)

	return acosh.(1 .+ (cosh(α * R) - 1) .* Ui) / α

end

"""
	function function sample_ϕ(N)

Function used to sample the angular coordinates ϕᵢ.
"""
function sample_ϕ(N::Number; distr = :uniform)
	return rand(Uniform(-π, π), N)
end

"""
	function hyper_d(r1, ϕ1, r2, ϕ2; ζ = 1)

Function used to calculate the hyperbolic distance between two nodes `i` and `j`:
- node `i`: `(r1, ϕ1)`
- node `j`: `(r2, ϕ2)`

Output:
- hyperbolic distance `x` between `(r1, ϕ1)` and `(r2, ϕ2)`.
"""
function hyper_d(r1, ϕ1, r2, ϕ2; ζ = 1)
	
	Δϕ = π - abs(π - abs(ϕ1 - ϕ2))
	k1 = cosh(ζ * r1) * cosh(ζ * r2)
	k2 = sinh(ζ * r1) * sinh(ζ * r2) * cos(Δϕ)
	
	return acosh(k1 - k2) / ζ
end

"""
	function connect_nodes_HRGG(r_i, ϕ_i, R; ζ = 1)

Function used to connect all nodes in a Hyperbolic RGG with T -> 0.
Output is the resulting adjacency matrix.

"""
function connect_nodes_HRGG(r_i, ϕ_i, R; ζ = 1)
	
	if length(r_i) != length(ϕ_i)
		error("Length mismatch r_i and ϕ_i")
	end
	
	N = length(r_i)
	A = zeros(N, N)
	
	for i in 1:N
		for j in i+1:N
			if hyper_d(r_i[i], ϕ_i[i], r_i[j], ϕ_i[j]; ζ = ζ) < R
				A[i,j] = 1
				A[j,i] = 1
			end
		end
	end

	return A
end

"""
	function connect_nodes_softHRGG(r_i, ϕ_i, R, T; ζ = 1)

Function used to connect all nodes in a soft Hyperbolic RGG.
Output is the resulting adjacency matrix.

For soft HRGG the connection probability is:

p(x) = 1 / (1 + exp(β * (ζ / 2) * (x - R)))

where `x` is the hyperbolic distance between `(r_1, ϕ_1)` and `(r_2, ϕ_2)`.

"""
function connect_nodes_softHRGG(r_i, ϕ_i, R, T; ζ = 1)
	
	if length(r_i) != length(ϕ_i)
		error("Length mismatch r_i and ϕ_i")
	end
	
	N = length(r_i)
	A = zeros(N, N)
	
	for i in 1:N
		for j in i+1:N
			
			x = hyper_d(r_i[i], ϕ_i[i], r_i[j], ϕ_i[j]; ζ = ζ)
			p = px_softHRGG(R, x, T, ζ)
			
			if rand() < p
				A[i,j] = 1
				A[j,i] = 1
			end
		end
	end

	return A
end