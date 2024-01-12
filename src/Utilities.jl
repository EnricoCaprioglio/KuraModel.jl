using LaTeXStrings
using Distributions
using Plots

"""
	function Laplacian(A::AbstractMatrix; method = :graph)

Function used to compute the Laplacian of the graph defined by the
adjacency matrix `A`.

If `method == :graph` --> compute the graph Laplacian;

If `method == :normalized` --> compute the normalized graph Laplacian.

"""
function Laplacian(A::AbstractMatrix; method = :graph)

	N = length(A[1, :])
	
	if method == :graph
		
		D = diagm([sum(A[i, :]) for i in 1:N])
		
		return D - A
		
	elseif method == :normalized

		𝐼 = diagm(ones(N))
		k = [sum(A[i, :]) for i in 1:N]
		A′ = zeros(N, N)
		
		for i in 1:N
			for j in 1:N
				if A[i,j] != 0
					A′[i,j] = A[i,j] / √(k[i]*k[j])
				end
			end
		end

		return 𝐼 - A′
	end
	
end

"""
    function get_splits(Nc)
This function is used in many other functions, such as ω_macro or θs_macro.
Given an Array containing the number of nodes in each community, such that
`Nc[i]` is the number of nodes in community `c_{i}` it outputs a useful iterator
to perform calculations for each community.

See functions `ω_macro` or `θs_macro` to see how it us used, or the example below
to find the average ω for each community.
## Example: 
```jldoctest
julia> Nc = [4,3,2,3]
splits = get_splits(Nc)
>> 8-element Vector{Any}:
  1
  4
  5
  7
  8
  9
 10
 12
 ω = [1,2,3,2,4,5,6,3,3,7,8,9]
 means=[]
 for i in 1:length(splits)
     if i%2 == 1
         push!(means,mean(ω[splits[i]:splits[i+1]]))
     end
 end
 means
>> 4-element Vector{Any}:
2.0
5.0
3.0
8.0
 ```
"""
function get_splits(Nc)
	
	splits = [1, Nc[1]]
	
	for i in 2:lastindex(Nc)
		
		first = sum(Nc[1:(i-1)]) + 1
		last = sum(Nc[1:i])
		
		push!(splits, first)
		push!(splits, last)
	end
	
	return splits
end
# function get_splits(Nc)
#     splits = []
# 	for (i,c) in enumerate(Nc)
# 		append!(splits,sum(Nc[1:i-1])+1, sum(Nc[1:i-1])+Nc[i])
# 	end
#     return splits
# end

"""
```jldoctest
	function plot_local_op(θs::AbstractArray, Nc::AbstractArray; steps_to_plot = nothing)
```
Function used to plot the local order parameters for each community `c`.
To plot use `plot(plot_local_op(θs, Nc)..., kwargs)`.
"""
function plot_local_op(θs::AbstractArray, Nc::AbstractArray; steps_to_plot = nothing)
	
    ps = []  # store plots
	splits = get_splits(Nc)

	for i in 1:lastindex(splits)
		if i % 2 == 1
			c = round(Integer, ceil(i/2))
			if steps_to_plot === nothing
				push!(ps,
					plot(
						macro_op(θs[ : , splits[i] : splits[i+1]]),
						label = "",
						ylims = [0., 1.],
						bg = RGB(0.2, 0.2, 0.2),
						title = L"R_{%$(c)}(t)"
					)
				)
			else
				push!(ps,
					plot(
						steps_to_plot,
						macro_op(θs[ : , splits[i] : splits[i+1]]),
						label = "",
						ylims = [0., 1.],
						bg = Plots.RGB(0.2, 0.2, 0.2),
						title = L"R_{%$(c)}(t)"
					)
				)
			end
		end
	end
	return ps
end

"""
    fun_pattern(A::AbstractMatrix, θs::AbstractArray, t_span)

Function used to calculate the functional pattern (correlation matrix)
of a system of oscillators.

    Inputs:
        A::AbstractMatrix   adjacency matrix of the underlying network
        θs::AbstractArray   instantaneous phases at each timestep
        t_span              time window of interest (can be a tuple or vector of order 2)

    Output:
        R::AbstractMatrix   functional pattern matrix for time window t_span
"""
function fun_pattern(A::AbstractMatrix, θs::AbstractArray, t_span)
    
    R = zeros(size(A))

    for i in 1:length(A[1,:])
        for j in 1:length(A[1,:])
            local corr_sum = 0
            for t in t_span[1]:t_span[2]
                corr_sum += cos(θs[t,j] - θs[t,i])
            end
            R[i,j] = corr_sum / (t_span[2] - t_span[1]+1)
        end
    end

    return R
end

"""
    function ω_macro(ω::AbstractArray, Nc::AbstractArray)
This function finds the local mean natural frequency for each community.

    Inputs:
        ω::AbstractArray		the vector of frequencies.
        Nc::AbstractArray       array containing the number of oscillators per community

    Output:
        means::AbstractVector   the natural frequencies of each macroscillator
"""
function ω_macro(ω::AbstractArray,Nc::AbstractArray)

    splits = get_splits(Nc)
	means = []
    
	for i in 1:length(splits)
        if i%2 == 1
            push!(means, mean(ω[splits[i]:splits[i+1]]))
        end
    end
	return means
end

"""
    function θs_macro(θs::AbstractArray,Nc::AbstractArray)
This function finds the local mean phase for each community at each time step.

    Inputs:
        θs::AbstractArray		    the matrix of instantaneous phases.
        Nc::AbstractArray           array containing the number of oscillators per community
    Output:
        θs_means::AbstractArray     instantaneous community averaged phases at each timestep

The output is a matrix (timesteps × M).
"""
function θs_macro(θs::AbstractArray,Nc::AbstractArray)
	
    tot_steps=length(θs[:,1])
	θs_means=zeros(tot_steps,length(Nc))
    splits = get_splits(Nc)

	for t in 1:tot_steps
        for i in 1:length(splits)
            if i%2 == 1
                c=round(Integer,ceil(i/2)) # find community index
                θs_means[t,c]=mean(θs[t,splits[i]:splits[i+1]])
            end
        end
	end

	return θs_means
end

"""
```jldoctest
	function ω_locals(ω::AbstractArray, Nc::AbstractArray)
```
this function finds the local mean natural frequency for each community.
"""
function ω_locals(ω::AbstractArray, Nc::AbstractArray)

	splits = get_splits(Nc)
	means=[]

	for i in 1:length(splits)
		if i%2 == 1
			push!(means, mean(ω[splits[i]:splits[i+1]]))
		end
	end

	return means
end

"""
```jldoctest
	function θs_locals(θs::AbstractArray, Nc::AbstractArray)
```
this function finds the local mean phase for each community at each time step.
The output is a matrix (timesteps × M)
"""
function θs_locals(θs::AbstractArray, Nc::AbstractArray)
	
	splits = get_splits(Nc)
	tot_steps = length(θs[:,1])
	θs_means = zeros(tot_steps, length(Nc))
	
	for t in 1:tot_steps
		for i in 1:length(splits)
			if i%2 == 1
				c_temp=round(Integer, ceil(i/2))
				θs_means[t,c_temp] = mean(θs[t, splits[i]:splits[i+1]])
			end
		end
	end
	return θs_means
end

"""
```jldoctest
	get_neighbours(
		A::AbstractMatrix,
		i::Integer;
		val = 0.0)
```
This function outputs the list of neighbours of node `i`.

The optional argument `val` is used to find the list of neighbours with a shared weighted edge of weight equal to `val`.

Parameters:

	A: 			adjacency matrix
	i: 			node index

Optional parameters:

	val: 		if the matrix is weighted

"""
function get_neighbours(A::AbstractMatrix, i::Integer; val = 0.0)
	
	if val != 0.0
		return findall(x -> x == val, A[i, :])
	else
		return findall(x -> x != 0.0, A[i, :])
	end
	
end

"""
```jldoctest
	macro_to_micro(
		Amacro::AbstractMatrix,
		C::Integer,
		M::Integer,
		d1::Integer;
		d0 = 0.0, μ = 0.6, ν = 0.4
	)
```
This function creates a community-structured netwrok with C communities of M nodes per community. Each node has `d0` intra-community connections and `d1` inter-community connections on average.

Parameters:

	Amacro: 	macroscopic connectivity matrix
	C: 			number of communities
	M: 			number of nodes per community
	d1: 		number of inter-community connections

Optional parameters:

	d0: 		number of intra-community connections
	μ: 			coupling between agents in the same community
	ν: 			coupling between agents in different communities

Above C × M = 2000 nodes the function will start to slow down.
"""
function macro_to_micro(Amacro::AbstractMatrix, C::Integer, M::Integer, d1::Real; d0 = 0.0, μ = 0.6, ν = 0.4)

	N = C * M
	# init matrix
	A = zeros(N, N)
	sp = get_splits(repeat([M], C))
		
	# add communities
	if d0 == 0.0
		for i in 1:lastindex(sp)
			if i % 2 == 1
				local mat = ones(M, M) .* μ
				A[sp[i] : sp[i+1], sp[i] : sp[i+1]] = mat
			end
		end
	else
		for i in 1:N
			# i belongs to community x1
	        x1 = mod(ceil(i / M) - 1, C) + 1
	        for j in (i + 1):N
	            if i != j
					# j belongs to community x2
	                x2 = mod(ceil(j / M) - 1, C) + 1
	                if x1 == x2 && rand() < (d0 / (M - 1))
	                    A[i, j] = μ
                    	A[j, i] = μ
					end
				end
			end
		end
	end
	
	# add inter-community connectivity
	for i in 1:N
		# i belongs to community x1
		x1 = round(Integer, mod(ceil(i / M) - 1, C) + 1)
		# neighbourhood of community x1:
		local macro_neigh = get_neighbours(Amacro, x1)
		for j in (i + 1):N
			# j belongs to community x2
			x2 = round(Integer, mod(ceil(j / M) - 1, C) + 1)
			if x1 != x2 && x2 ∈ macro_neigh  # j is a possible connection for i
				if rand() < (d1 / (M * length(macro_neigh)))
					A[i, j] = ν
                    A[j, i] = ν
				end
			end
		end	
	end

	setindex!.(Ref(A), 0.0, 1:N, 1:N)  # set diagonal to zero
	
	return A
end